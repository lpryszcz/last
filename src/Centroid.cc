// Copyright 2008, 2009, 2010, 2011 Michiaki Hamada
// Copyright 2012, 2013 Toshiyuki Sato

#include "Centroid.hh"
#include "GappedXdropAlignerInl.hh"
#include <algorithm>
#include <cassert>
#include <cmath> // for exp
#include <cfloat>   // for DBL_MAX
#include <cstdlib>  // for abs
#include <iomanip>

#define CI(type) std::vector<type>::const_iterator  // added by MCF

static const double DINF = DBL_MAX / 2;

namespace{
  double EXP ( double x ) {
    return std::exp (x);
  }
}


namespace cbrc{

  ExpectedCount::ExpectedCount ()
  {
    double d0 = 0;
    MM = d0; MD = d0; MP = d0; MI = d0; MQ = d0;
    DD = d0; DM = d0; DI = d0;
    PP = d0; PM = d0; PD = d0; PI = d0;
    II = d0; IM = d0;
    SM = d0; SD = d0; SP = d0; SI = d0; SQ = d0;

    for (int n=0; n<scoreMatrixRowSize; n++)
      for (int m=0; m<scoreMatrixRowSize; m++) emit[n][m] = d0;
  }

  std::ostream& ExpectedCount::write (std::ostream& os) const
  {
    for (int n=0; n<scoreMatrixRowSize; ++n) {
      for (int m=0; m<scoreMatrixRowSize; ++m) {
	double prob = emit[n][m];
	if (prob > 0)
	  os << "emit[" << n << "][" << m << "]=" << emit[n][m] << std::endl;
      }
    }
    os << "M->M=" << MM << std::endl;
    os << "M->D=" << MD << std::endl;
    os << "M->P=" << MP << std::endl;
    os << "M->I=" << MI << std::endl;

    os << "D->D=" << DD << std::endl;
    os << "D->M=" << DM << std::endl;
    os << "D->I=" << DI << std::endl;

    os << "P->P=" << PP << std::endl;
    os << "P->M=" << PM << std::endl;
    os << "P->D=" << PD << std::endl;
    os << "P->I=" << PI << std::endl;

    os << "I->I=" << II << std::endl;
    os << "I->M=" << IM << std::endl;

    return os;
  }

  void Centroid::setScoreMatrix( const ScoreMatrixRow* sm, double T ) {
    this -> T = T;
    this -> isPssm = false;
    for ( int n=0; n<scoreMatrixRowSize; ++n )
      for ( int m=0; m<scoreMatrixRowSize; ++m ) {
	match_score[n][m] = EXP ( sm[ n ][ m ] / T );
      }
  }

  void Centroid::setPssm( const ScoreMatrixRow* pssm, size_t qsize, double T,
			  const OneQualityExpMatrix& oqem,
			  const uchar* sequenceBeg, const uchar* qualityBeg ) {
    this->T = T;
    this -> isPssm = true;
    pssmExp.resize( qsize * scoreMatrixRowSize );
    pssmExp2 = reinterpret_cast<ExpMatrixRow*> ( &pssmExp[0] );

    if( oqem ){  // fast special case
      makePositionSpecificExpMatrix( oqem, sequenceBeg, sequenceBeg + qsize,
                                     qualityBeg, &pssmExp[0] );
    }
    else{  // slow general case
      for ( size_t i=0; i<qsize; ++i ) {
        for ( unsigned j=0; j<scoreMatrixRowSize; ++j ) {
          pssmExp2[ i ][ j ] = EXP ( pssm[ i ][ j ] / T );
        }
      }
    }
  }

  void Centroid::initForwardMatrix(){
    scale.assign ( numAntidiagonals + 2, 1.0 ); // scaling
    size_t n = xa.scoreEndIndex( numAntidiagonals );

    if ( fM.size() < n ) {
      fM.resize( n );
      fD.resize( n );
      fI.resize( n );
      fP.resize( n );
    }

    fM[0] = 1;
  }

  void Centroid::initBackwardMatrix(){
    mD.assign( numAntidiagonals + 2, 0.0 );
    mI.assign( numAntidiagonals + 2, 0.0 );

    size_t n = xa.scoreEndIndex( numAntidiagonals );
    bM.assign( n, 0.0 );
    bD.assign( n, 0.0 );
    bI.assign( n, 0.0 );
    bP.assign( n, 0.0 );
  }

  void Centroid::initDecodingMatrix(){
    X.resize( fM.size() );
  }

  void Centroid::updateScore( double score, size_t antiDiagonal, size_t cur ){
    if( bestScore < score ){
      bestScore = score;
      bestAntiDiagonal = antiDiagonal;
      bestPos1 = cur;
    }
  }

  // xxx this will go wrong for non-delimiters with severe mismatch scores
  static bool isDelimiter(uchar c, const double *expScores) {
    return expScores[c] <= 0.0;
  }

  void Centroid::forward( const uchar* seq1, const uchar* seq2,
			  size_t start1, size_t start2,
			  bool isForward, int globality,
			  const GeneralizedAffineGapCosts& gap ){

    //std::cout << "[forward] start1=" << start1 << "," << "start2=" << start2 << "," << "isForward=" << isForward << std::endl;
    seq1 += start1;
    seq2 += start2;
    const ExpMatrixRow* pssm = isPssm ? pssmExp2 + start2 : 0;
    const int seqIncrement = isForward ? 1 : -1;

    initForwardMatrix();

    double Z = 0.0;  // partion function of forward values

    const bool isAffine = gap.isAffine();
    const int E = gap.delExtend;
    const int F = gap.delExist;
    const int EI = gap.insExtend;
    const int FI = gap.insExist;
    const int P = gap.pairExtend;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eEI = EXP ( - EI / T );
    const double eFI = EXP ( - FI / T );
    const double eP = EXP ( - P / T );

    assert( gap.insExist == gap.delExist || eP <= 0.0 );

    for( size_t k = 0; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t seq1beg = seq1start( k );
      const size_t seq2pos = k - seq1beg;
      const double scale12 = scale[k+1] * scale[k];
      const double scale1  = scale[k+1];
      double sum_f = 0.0; // sum of forward values

      const double seE = eE * scale1;
      const double seEI = eEI * scale1;
      const double seP = eP * scale12;

      const size_t scoreEnd = xa.scoreEndIndex( k );
      double* fM0 = &fM[ scoreEnd ];
      double* fD0 = &fD[ scoreEnd ];
      double* fI0 = &fI[ scoreEnd ];
      double* fP0 = &fP[ scoreEnd ];

      const size_t horiBeg = xa.hori( k, seq1beg );
      const size_t vertBeg = xa.vert( k, seq1beg );
      const size_t diagBeg = xa.diag( k, seq1beg );
      const double* fD1 = &fD[ horiBeg ];
      const double* fI1 = &fI[ vertBeg ];
      const double* fM2 = &fM[ diagBeg ];
      const double* fP2 = &fP[ diagBeg ];

      const double* fM0last = fM0 + xa.numCellsAndPads( k ) - 1;

      const uchar* s1 = seqPtr( seq1, isForward, seq1beg );

      *fM0++ = *fD0++ = *fI0++ = *fP0++ = 0.0;  // add one pad cell

      if (! isPssm) {
	const uchar* s2 = seqPtr( seq2, isForward, seq2pos );

	if (isAffine) {
	  while (1) {
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    *fD0 = xM * eF + xD;
	    *fI0 = (xM + xD) * eFI + xI;
	    *fM0 = (xM + xD + xI) * match_score[ *s1 ][ *s2 ];
	    sum_f += xM;
	    if( globality && (isDelimiter(*s2, *match_score) ||
			      isDelimiter(*s1, *match_score)) ){
	      Z += xM + xD + xI;
	    }
	    if (fM0 == fM0last) break;
	    fM0++; fD0++; fI0++;
	    fM2++; fD1++; fI1++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	  }
	} else {
	  while (1) {
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    const double xP = *fP2 * seP;
	    *fD0 = xM * eF + xD + xP;
	    *fI0 = (xM + xD) * eFI + xI + xP;
	    *fM0 = (xM + xD + xI + xP) * match_score[ *s1 ][ *s2 ];
	    *fP0 = xM * eF + xP;
	    sum_f += xM;
	    if( globality && (isDelimiter(*s2, *match_score) ||
			      isDelimiter(*s1, *match_score)) ){
	      Z += xM + xD + xI + xP;
	    }
	    if (fM0 == fM0last) break;
	    fM0++; fD0++; fI0++; fP0++;
	    fM2++; fD1++; fI1++; fP2++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	  }
	}
      } else {
	const ExpMatrixRow* p2 = seqPtr( pssm, isForward, seq2pos );

	if (isAffine) {
	  while (1) {
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    *fD0 = xM * eF + xD;
	    *fI0 = (xM + xD) * eFI + xI;
	    *fM0 = (xM + xD + xI) * (*p2)[ *s1 ];
	    sum_f += xM;
	    if ( globality && (isDelimiter(0, *p2) ||
			       isDelimiter(*s1, *pssm)) ){
	      Z += xM + xD + xI;
	    }
	    if (fM0 == fM0last) break;
	    fM0++; fD0++; fI0++;
	    fM2++; fD1++; fI1++;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }
	} else {
	  while (1) {
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    const double xP = *fP2 * seP;
	    *fD0 = xM * eF + xD + xP;
	    *fI0 = (xM + xD) * eFI + xI + xP;
	    *fM0 = (xM + xD + xI + xP) * (*p2)[ *s1 ];
	    *fP0 = xM * eF + xP;
	    sum_f += xM;
	    if ( globality && (isDelimiter(0, *p2) ||
			       isDelimiter(*s1, *pssm)) ){
	      Z += xM + xD + xI + xP;
	    }
	    if (fM0 == fM0last) break;
	    fM0++; fD0++; fI0++; fP0++;
	    fM2++; fD1++; fI1++; fP2++;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }
	}
      }

      if( !globality ) Z += sum_f;
      scale[k+2] = 1.0 / (sum_f + 1.0);  // seems ugly
      Z *= scale[k+2]; // scaling
    } // k

    //std::cout << "# Z=" << Z << std::endl;
    assert( Z > 0.0 );
    scale[ numAntidiagonals + 1 ] /= Z;  // this causes scaled Z to equal 1
  }

  // added by M. Hamada
  // compute posterior probabilities while executing backward algorithm
  void Centroid::backward( const uchar* seq1, const uchar* seq2,
			   size_t start1, size_t start2,
			   bool isForward, int globality,
			   const GeneralizedAffineGapCosts& gap ){

    //std::cout << "[backward] start1=" << start1 << "," << "start2=" << start2 << "," << "isForward=" << isForward << std::endl;
    seq1 += start1;
    seq2 += start2;
    const ExpMatrixRow* pssm = isPssm ? pssmExp2 + start2 : 0;
    const int seqIncrement = isForward ? 1 : -1;

    initBackwardMatrix();

    const bool isAffine = gap.isAffine();
    const int E = gap.delExtend;
    const int F = gap.delExist;
    const int EI = gap.insExtend;
    const int FI = gap.insExist;
    const int P = gap.pairExtend;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eEI = EXP ( - EI / T );
    const double eFI = EXP ( - FI / T );
    const double eP = EXP ( - P / T );
    double scaledUnit = 1.0;

    assert( gap.insExist == gap.delExist || eP <= 0.0 );

    for( size_t k = numAntidiagonals; k-- > 0; ){
      const size_t seq1beg = seq1start( k );
      const size_t seq2pos = k - seq1beg;
      const double scale12 = scale[k+1] * scale[k];
      const double scale1  = scale[k+1];
      scaledUnit *= scale[k+2];

      const double seE = eE * scale1;
      const double seEI = eEI * scale1;
      const double seP = eP * scale12;

      const size_t scoreEnd = xa.scoreEndIndex( k );
      const double* bM0 = &bM[ scoreEnd + 1 ];
      const double* bD0 = &bD[ scoreEnd + 1 ];
      const double* bI0 = &bI[ scoreEnd + 1 ];
      const double* bP0 = &bP[ scoreEnd + 1 ];

      const size_t horiBeg = xa.hori( k, seq1beg );
      const size_t vertBeg = xa.vert( k, seq1beg );
      const size_t diagBeg = xa.diag( k, seq1beg );
      double* bD1 = &bD[ horiBeg ];
      double* bI1 = &bI[ vertBeg ];
      double* bM2 = &bM[ diagBeg ];
      double* bP2 = &bP[ diagBeg ];

      const double* fD1 = &fD[ horiBeg ];
      const double* fI1 = &fI[ vertBeg ];
      const double* fP2 = &fP[ diagBeg ];

      const double* bM0last = bM0 + xa.numCellsAndPads( k ) - 2;

      double* mDout = &mD[ seq1beg ];
      double* mIout = &mI[ seq2pos ];

      const uchar* s1 = seqPtr( seq1, isForward, seq1beg );

      if (! isPssm ) {
	const uchar* s2 = seqPtr( seq2, isForward, seq2pos );

	if (isAffine) {
	  while (1) {
	    double yM = *bM0 * match_score[ *s1 ][ *s2 ];
	    double yD = *bD0;
	    double yI = *bI0;
	    double zM = yM + yD * eF + yI * eFI;
	    double zD = yM + yD + yI * eFI;
	    double zI = yM + yI;
	    if( globality ){
	      if( isDelimiter(*s2, *match_score) ||
		  isDelimiter(*s1, *match_score) ){
		zM += scaledUnit;  zD += scaledUnit;  zI += scaledUnit;
	      }
	    }else{
	      zM += scaledUnit;
	    }
	    *bM2 = zM * scale12;
	    *bD1 = zD * seE;
	    *bI1 = zI * seEI;

	    *mDout += (*fD1) * (*bD1);
	    *mIout += (*fI1) * (*bI1);

	    if (bM0 == bM0last) break;
	    mDout++; mIout--;
	    bM2++; bD1++; bI1++;
	    bM0++; bD0++; bI0++;
	    fD1++; fI1++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	  }
	} else {
	  while (1) {
	    double yM = *bM0 * match_score[ *s1 ][ *s2 ];
	    double yD = *bD0;
	    double yI = *bI0;
	    double yP = *bP0;
	    double zM = yM + yD * eF + yI * eFI + yP * eF;
	    double zD = yM + yD + yI * eFI;
	    double zI = yM + yI;
	    double zP = yM + yP + yD + yI;
	    if( globality ){
	      if( isDelimiter(*s2, *match_score) ||
		  isDelimiter(*s1, *match_score) ){
		zM += scaledUnit;  zD += scaledUnit;  zI += scaledUnit;
		zP += scaledUnit;
	      }
	    }else{
	      zM += scaledUnit;
	    }
	    *bM2 = zM * scale12;
	    *bD1 = zD * seE;
	    *bI1 = zI * seEI;
	    *bP2 = zP * seP;

	    double probp = *fP2 * *bP2;
	    *mDout += (*fD1) * (*bD1) + probp;
	    *mIout += (*fI1) * (*bI1) + probp;

	    if (bM0 == bM0last) break;
	    mDout++; mIout--;
	    bM2++; bD1++; bI1++; bP2++;
	    bM0++; bD0++; bI0++; bP0++;
	    fD1++; fI1++; fP2++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	  }
	}
      } else {
	const ExpMatrixRow* p2 = seqPtr( pssm, isForward, seq2pos );

	if (isAffine) {
	  while (1) {
	    double yM = *bM0 * ( *p2 )[ *s1 ];
	    double yD = *bD0;
	    double yI = *bI0;
	    double zM = yM + yD * eF + yI * eFI;
	    double zD = yM + yD + yI * eFI;
	    double zI = yM + yI;
	    if( globality ){
	      if( isDelimiter(0, *p2) ||
		  isDelimiter(*s1, *pssm) ){
		zM += scaledUnit;  zD += scaledUnit;  zI += scaledUnit;
	      }
	    }else{
	      zM += scaledUnit;
	    }
	    *bM2 = zM * scale12;
	    *bD1 = zD * seE;
	    *bI1 = zI * seEI;

	    *mDout += (*fD1) * (*bD1);
	    *mIout += (*fI1) * (*bI1);

	    if (bM0 == bM0last) break;
	    mDout++; mIout--;
	    bM2++; bD1++; bI1++;
	    bM0++; bD0++; bI0++;
	    fD1++; fI1++;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }
	} else {
	  while (1) {
	    double yM = *bM0 * ( *p2 )[ *s1 ];
	    double yD = *bD0;
	    double yI = *bI0;
	    double yP = *bP0;
	    double zM = yM + yD * eF + yI * eFI + yP * eF;
	    double zD = yM + yD + yI * eFI;
	    double zI = yM + yI;
	    double zP = yM + yP + yD + yI;
	    if( globality ){
	      if( isDelimiter(0, *p2) ||
		  isDelimiter(*s1, *pssm) ){
		zM += scaledUnit;  zD += scaledUnit;  zI += scaledUnit;
		zP += scaledUnit;
	      }
	    }else{
	      zM += scaledUnit;
	    }
	    *bM2 = zM * scale12;
	    *bD1 = zD * seE;
	    *bI1 = zI * seEI;
	    *bP2 = zP * seP;

	    double probp = *fP2 * *bP2;
	    *mDout += (*fD1) * (*bD1) + probp;
	    *mIout += (*fI1) * (*bI1) + probp;

	    if (bM0 == bM0last) break;
	    mDout++; mIout--;
	    bM2++; bD1++; bI1++; bP2++;
	    bM0++; bD0++; bI0++; bP0++;
	    fD1++; fI1++; fP2++;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }
	}
      }
    }

    //std::cout << "# bM[0]=" << bM[0] << std::endl;
    //ExpectedCount ec;
    //computeExpectedCounts ( seq1, seq2, start1, start2, isForward, gap, ec );
    //ec.write (std::cerr);
  }

  double Centroid::dp( double gamma ){
    if (outputType == 5 ) return dp_centroid( gamma );
    else if (outputType == 6 ) return dp_ama (gamma);
    return 0;
  }

  void Centroid::traceback( std::vector< SegmentPair >& chunks,
			    double gamma ) const{
    if (outputType==5) traceback_centroid( chunks, gamma );
    else if (outputType==6) traceback_ama( chunks, gamma);
  }

  double Centroid::dp_centroid( double gamma ){

    initDecodingMatrix();

    for( size_t k = 1; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      size_t seq1pos = seq1start( k );

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const size_t h = xa.hori( k, seq1pos );
      const size_t d = xa.diag( k, seq1pos );
      const double* X1 = &X[h];
      const double* X2 = &X[d];
      const double* fM2 = &fM[d];
      const double* bM2 = &bM[d];

      *X0++ = -DINF;		// add one pad cell

      do{
	const double matchProb = (*fM2++) * (*bM2++);
	const double s = ( gamma + 1 ) * matchProb - 1;
	const double oldX1 = *X1++;  // Added by MCF
	const double score = std::max( std::max( oldX1, *X1 ), *X2++ + s );
	//assert ( score >= 0 );
	updateScore ( score, k, seq1pos );
	*X0++ = score;
	seq1pos++;
      }while( X0 != x0end );
    }
    return bestScore;
  }

  void Centroid::traceback_centroid( std::vector< SegmentPair >& chunks,
				     double gamma ) const{
    //std::cout << "[c] bestAntiDiagonal=" << bestAntiDiagonal << ": bestPos1=" << bestPos1 << std::endl;

    size_t k = bestAntiDiagonal;
    size_t i = bestPos1;
    size_t oldPos1 = i;

    while( k > 0 ){
      const size_t h = xa.hori( k, i );
      const size_t v = xa.vert( k, i );
      const size_t d = xa.diag( k, i );
      const double matchProb = fM[d] * bM[d];
      const int m = maxIndex( X[d] + (gamma + 1) * matchProb - 1, X[h], X[v] );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 0 ){
	chunks.push_back( SegmentPair( i, k - i, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  double Centroid::dp_ama( double gamma ){

    initDecodingMatrix();

    mX1.assign ( numAntidiagonals + 2, 1.0 );
    mX2.assign ( numAntidiagonals + 2, 1.0 );

    for (size_t k = 0; k < numAntidiagonals; ++k) {
      size_t seq1pos = seq1start(k);
      size_t seq2pos = k - seq1pos;
      size_t loopBeg = xa.diag(k, seq1pos);
      size_t loopEnd = loopBeg + xa.numCellsAndPads(k) - 1;
      for (size_t i = loopBeg; i < loopEnd; ++i) {
	const double matchProb = fM[i] * bM[i];
	mX1[seq1pos++] -= matchProb;
	mX2[seq2pos--] -= matchProb;
      }
    }

    for( size_t k = 1; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      size_t seq1pos = seq1start( k );
      size_t seq2pos = k - seq1pos;

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const size_t h = xa.hori( k, seq1pos );
      const size_t d = xa.diag( k, seq1pos );
      const double* X1 = &X[h];
      const double* X2 = &X[d];
      const double* fM2 = &fM[d];
      const double* bM2 = &bM[d];

      *X0++ = -DINF;		// add one pad cell

      do{
	const double matchProb = (*fM2++) * (*bM2++);
	const double thisD = mD[seq1pos];
	const double thisI = mI[seq2pos];
	const double thisXD = mX1[seq1pos] - thisD;
	const double thisXI = mX2[seq2pos] - thisI;
	const double s = 2 * gamma * matchProb - (thisXD + thisXI);
	const double u = gamma * thisD - thisXD;
	const double t = gamma * thisI - thisXI;
	const double oldX1 = *X1++;  // Added by MCF
	const double score = std::max( std::max( oldX1 + u, *X1 + t), *X2++ + s );
	updateScore ( score, k, seq1pos );
	*X0++ = score;
	seq1pos++;
	seq2pos--;
      }while( X0 != x0end );
    }

    return bestScore;
  }

  void Centroid::traceback_ama( std::vector< SegmentPair >& chunks,
			    double gamma ) const{
    //std::cout << "[c] bestAntiDiagonal=" << bestAntiDiagonal << ": bestPos1=" << bestPos1 << std::endl;

    size_t k = bestAntiDiagonal;
    size_t i = bestPos1;
    size_t oldPos1 = i;

    while( k > 0 ){
      const size_t j = k - i;
      const size_t h = xa.hori( k, i );
      const size_t v = xa.vert( k, i );
      const size_t d = xa.diag( k, i );
      const double matchProb = fM[d] * bM[d];
      const double thisD = mD[i];
      const double thisI = mI[j];
      const double thisXD = mX1[i] - thisD;
      const double thisXI = mX2[j] - thisI;
      const double s = 2 * gamma * matchProb - (thisXD + thisXI);
      const double t = gamma * thisI - thisXI;
      const double u = gamma * thisD - thisXD;
      const int m = maxIndex( X[d] + s, X[h] + u, X[v] + t );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 0 ){
	chunks.push_back( SegmentPair( i, k - i, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  // Return an ASCII code representing an error probability.  The
  // printable codes are 33--126, but here we use a maximum of 125, so
  // that 126 is reserved for special cases.
  static uchar asciiProbability( double probCorrect ){
    assert( probCorrect >= 0 );
    //assert( probCorrect <= 1 );  // can fail: floating point is imperfect
    double e = 1 - probCorrect;
    double f = std::max( e, 1e-10 );  // avoid overflow errors
    double g = -10 * std::log10(f);
    int i = static_cast<int>(g);  // round fractions down
    int j = i + 33;
    int k = std::min( j, 125 );
    return static_cast<uchar>(k);
  }

  void Centroid::getMatchAmbiguities(std::vector<uchar>& ambiguityCodes,
				     size_t seq1end, size_t seq2end,
				     size_t length) const {
    while (length) {
      size_t d = xa.diag(seq1end + seq2end, seq1end);
      double p = fM[d] * bM[d];
      ambiguityCodes.push_back(asciiProbability(p));
      --seq1end;  --seq2end;  --length;
    }
  }

  void Centroid::getDeleteAmbiguities(std::vector<uchar>& ambiguityCodes,
				      size_t seq1end, size_t seq1beg) const {
    for (size_t i = seq1end; i > seq1beg; --i)
      ambiguityCodes.push_back(asciiProbability(mD[i]));
  }

  void Centroid::getInsertAmbiguities(std::vector<uchar>& ambiguityCodes,
				      size_t seq2end, size_t seq2beg) const {
    for (size_t i = seq2end; i > seq2beg; --i)
      ambiguityCodes.push_back(asciiProbability(mI[i]));
  }

  double Centroid::logPartitionFunction() const{
    double x = 0.0;
    for( size_t k = 0; k < numAntidiagonals; ++k ){
      x -= std::log( scale[k+2] );
    }
    return T * x;
  }

  void Centroid::computeExpectedCounts ( const uchar* seq1, const uchar* seq2,
					 size_t start1, size_t start2,
					 bool isForward,
					 const GeneralizedAffineGapCosts& gap,
					 ExpectedCount& c ) const{
    seq1 += start1;
    seq2 += start2;
    const ExpMatrixRow* pssm = isPssm ? pssmExp2 + start2 : 0;
    const int seqIncrement = isForward ? 1 : -1;

    const bool isAffine = gap.isAffine();
    const int E = gap.delExtend;
    const int F = gap.delExist;
    const int EI = gap.insExtend;
    const int FI = gap.insExist;
    const int P = gap.pairExtend;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eEI = EXP ( - EI / T );
    const double eFI = EXP ( - FI / T );
    const double eP = EXP ( - P / T );

    assert( gap.insExist == gap.delExist || eP <= 0.0 );

    for( size_t k = 0; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t seq1beg = seq1start( k );
      const size_t seq2pos = k - seq1beg;
      const double scale12 = scale[k+1] * scale[k];
      const double scale1  = scale[k+1];

      const double seE = eE * scale1;
      const double seEI = eEI * scale1;
      const double seP = eP * scale12;

      const uchar* s1 = seqPtr( seq1, isForward, seq1beg );
      const uchar* s2 = seqPtr( seq2, isForward, seq2pos );

      const size_t scoreEnd = xa.scoreEndIndex( k );
      const double* bM0 = &bM[ scoreEnd + 1 ];
      const double* bD0 = &bD[ scoreEnd + 1 ];
      const double* bI0 = &bI[ scoreEnd + 1 ];
      const double* bP0 = &bP[ scoreEnd + 1 ];

      const size_t horiBeg = xa.hori( k, seq1beg );
      const size_t vertBeg = xa.vert( k, seq1beg );
      const size_t diagBeg = xa.diag( k, seq1beg );
      const double* fD1 = &fD[ horiBeg ];
      const double* fI1 = &fI[ vertBeg ];
      const double* fM2 = &fM[ diagBeg ];
      const double* fP2 = &fP[ diagBeg ];

      const double* bM0last = bM0 + xa.numCellsAndPads( k ) - 2;

      if (! isPssm ) {
	if (isAffine) {
	  while (1) {
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    const double yM = *bM0 * match_score[ *s1 ][ *s2 ];
	    const double yD = *bD0;
	    const double yI = *bI0;
	    c.emit[*s1][*s2] += (xM + xD + xI) * yM;
	    c.MM += xM * yM;
	    c.DM += xD * yM;
	    c.IM += xI * yM;
	    c.MD += xM * yD * eF;
	    c.DD += xD * yD;
	    c.MI += xM * yI * eFI;
	    c.DI += xD * yI * eFI;
	    c.II += xI * yI;
	    if (bM0 == bM0last) break;
	    fM2++; fD1++; fI1++;
	    bM0++; bD0++; bI0++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	  }
	} else {
	  while (1) {
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    const double xP = *fP2 * seP;
	    const double yM = *bM0 * match_score[ *s1 ][ *s2 ];
	    const double yD = *bD0;
	    const double yI = *bI0;
	    const double yP = *bP0;
	    c.emit[*s1][*s2] += (xM + xD + xI + xP) * yM;
	    c.MM += xM * yM;
	    c.DM += xD * yM;
	    c.IM += xI * yM;
	    c.PM += xP * yM;
	    c.MD += xM * yD * eF;
	    c.DD += xD * yD;
	    c.PD += xP * yD;
	    c.MI += xM * yI * eFI;
	    c.DI += xD * yI * eFI;
	    c.II += xI * yI;
	    c.PI += xP * yI;
	    c.MP += xM * yP * eF;
	    c.PP += xP * yP;
	    if (bM0 == bM0last) break;
	    fM2++; fD1++; fI1++; fP2++;
	    bM0++; bD0++; bI0++; bP0++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	  }
	}
      }
      else {
	const ExpMatrixRow* p2 = seqPtr( pssm, isForward, seq2pos );

	if (isAffine) {
	  while (1) { // inner most loop
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    const double yM = *bM0 * ( *p2 )[ *s1 ];
	    const double yD = *bD0;
	    const double yI = *bI0;
	    c.emit[*s1][*s2] += (xM + xD + xI) * yM;  // xxx
	    c.MM += xM * yM;
	    c.DM += xD * yM;
	    c.IM += xI * yM;
	    c.MD += xM * yD * eF;
	    c.DD += xD * yD;
	    c.MI += xM * yI * eFI;
	    c.DI += xD * yI * eFI;
	    c.II += xI * yI;
	    if (bM0 == bM0last) break;
	    fM2++; fD1++; fI1++;
	    bM0++; bD0++; bI0++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	    p2 -= seqIncrement;
	  }
	}else{
	  while (1) { // inner most loop
	    const double xM = *fM2 * scale12;
	    const double xD = *fD1 * seE;
	    const double xI = *fI1 * seEI;
	    const double xP = *fP2 * seP;
	    const double yM = *bM0 * ( *p2 )[ *s1 ];
	    const double yD = *bD0;
	    const double yI = *bI0;
	    const double yP = *bP0;
	    c.emit[*s1][*s2] += (xM + xD + xI + xP) * yM;  // xxx
	    c.MM += xM * yM;
	    c.DM += xD * yM;
	    c.IM += xI * yM;
	    c.PM += xP * yM;
	    c.MD += xM * yD * eF;
	    c.DD += xD * yD;
	    c.PD += xP * yD;
	    c.MI += xM * yI * eFI;
	    c.DI += xD * yI * eFI;
	    c.II += xI * yI;
	    c.PI += xP * yI;
	    c.MP += xM * yP * eF;
	    c.PP += xP * yP;
	    if (bM0 == bM0last) break;
	    fM2++; fD1++; fI1++; fP2++;
	    bM0++; bD0++; bI0++; bP0++;
	    s1 += seqIncrement;
	    s2 -= seqIncrement;
	    p2 -= seqIncrement;
	  }
	}
      }
    }
  }
}  // end namespace cbrc
