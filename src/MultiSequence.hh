// Copyright 2008, 2009, 2010, 2011, 2014 Martin C. Frith

// This class holds multiple sequences and their names.  The sequences
// are concatenated, with delimiters between them.

// The final sequence may be "unfinished".  This happens if the
// sequence data hits a memory limit before we finish reading it.

#ifndef MULTISEQUENCE_HH
#define MULTISEQUENCE_HH

#include "ScoreMatrixRow.hh"
#include "VectorOrMmap.hh"

#include <string>
#include <iosfwd>

namespace cbrc{

class MultiSequence{
 public:
  typedef LAST_INT_TYPE indexT;
  typedef unsigned char uchar;

  // initialize with leftmost delimiter pad, ready for appending sequences
  void initForAppending( indexT padSizeIn );

  // re-initialize, but keep the last sequence if it is unfinished
  void reinitForAppending();

  // read seqCount finished sequences, and their names, from binary files
  void fromFiles( const std::string& baseName, indexT seqCount,
                  size_t qualitiesPerLetter );

  // write all the finished sequences and their names to binary files
  void toFiles( const std::string& baseName ) const;

  // Append a sequence with delimiters.  Don't let the total size of
  // the concatenated sequences plus pads exceed maxSeqLen: thus it
  // may not finish reading the sequence.
  std::istream& appendFromFasta( std::istream& stream, indexT maxSeqLen );

  // As above, but read quality scores too.
  std::istream& appendFromFastq( std::istream& stream, indexT maxSeqLen );

  // As above, but read quality scores too.
  std::istream& appendFromPrb( std::istream& stream, indexT maxSeqLen,
			       unsigned alphSize, const uchar decode[] );

  // As above, but read a PSSM too, in PSI-BLAST ASCII format.
  std::istream& appendFromPssm( std::istream& stream, indexT maxSeqLen,
                                const uchar* lettersToNumbers,
                                bool isMaskLowercase );

  // finish the last sequence: add final pad and end coordinate
  void finish();

  // unfinish the last sequence: remove final pad and end coordinate
  void unfinish();

  // did we finish reading the last sequence?
  bool isFinished() const{ return ends.size() == nameEnds.size(); }

  // how many finished sequences are there?
  indexT finishedSequences() const{ return ends.size() - 1; }

  // total length of finished sequences plus delimiters
  indexT finishedSize() const{ return ends.back(); }

  // total length of finished and unfinished sequences plus delimiters
  size_t unfinishedSize() const{ return seq.size(); }

  // which sequence is the coordinate in?
  indexT whichSequence( indexT coordinate ) const;

  indexT padBeg( indexT seqNum ) const{ return ends[seqNum] - padSize; }
  indexT seqBeg( indexT seqNum ) const{ return ends[seqNum]; }
  indexT seqEnd( indexT seqNum ) const{ return ends[seqNum+1] - padSize; }
  indexT padEnd( indexT seqNum ) const{ return ends[seqNum+1]; }
  indexT seqLen( indexT seqNum ) const{ return seqEnd(seqNum)-seqBeg(seqNum); }
  indexT padLen( indexT seqNum ) const{ return padEnd(seqNum)-padBeg(seqNum); }
  std::string seqName( indexT seqNum ) const;

  // get a pointer to the start of the sequence data
  const uchar* seqReader() const{ return seq.begin(); }
  /***/ uchar* seqWriter()      { return &seq.v[0]; }

  // get a pointer to the start of the PSSM, or NULL if there is no PSSM
  // I am not totally sure about the reinterpret_cast...

  const ScoreMatrixRow* pssmReader() const{
    return pssm.empty() ? 0
        : reinterpret_cast< const ScoreMatrixRow* >( &pssm[0] );
  }

  /* */ ScoreMatrixRow* pssmWriter()      {
    return pssm.empty() ? 0
        : reinterpret_cast<       ScoreMatrixRow* >( &pssm[0] );
  }

  // get a pointer to the start of the quality data
  const uchar* qualityReader() const{ return qualityScores.begin(); }
  /***/ uchar* qualityWriter()      { return &qualityScores.v[0]; }

  // How many quality scores are there per letter?  There might be
  // none at all, one per letter, or several (e.g. 4) per letter.
  size_t qualsPerLetter() const { return qualityScoresPerLetter; }

 private:
  indexT padSize;  // number of delimiter chars between sequences
  VectorOrMmap<uchar> seq;  // concatenated sequences
  VectorOrMmap<indexT> ends;  // coordinates of ends of delimiter pads
  VectorOrMmap<char> names;  // concatenated sequence names (to save memory)
  VectorOrMmap<indexT> nameEnds;  // endpoints of the names

  std::vector<int> pssm;  // position-specific scoring matrix
  std::vector<uchar> pssmColumnLetters;  // which input column is which letter

  // The quality scores may be ASCII-coded: to get the real scores,
  // subtract e.g. 33 or 64.  The real scores might be related to
  // error probabilities in one of these ways:
  // Qphred = -10*log10(p)
  // Qsolexa = -10*log10(p/(1-p))
  VectorOrMmap<uchar> qualityScores;
  size_t qualityScoresPerLetter;

  // read a FASTA header: read the whole line but store just the first word
  std::istream& readFastaName( std::istream& stream );

  // read the letters above PSSM columns, so we know which column is which
  std::istream& readPssmHeader( std::istream& stream );

  // add a new sequence name
  void addName( std::string& name );

  // can we finish the last sequence and stay within the memory limit?
  bool isFinishable( indexT maxSeqLen ) const;
};

// Divide the sequences into a given number of roughly-equally-sized
// chunks, and return the first sequence in the Nth chunk.
inline size_t firstSequenceInChunk(const MultiSequence &m,
				   size_t numOfChunks, size_t chunkNum) {
  size_t numOfSeqs = m.finishedSequences();
  size_t beg = m.seqBeg(0);
  size_t end = m.padEnd(numOfSeqs - 1) - 1;
  unsigned long long len = end - beg;  // try to avoid overflow
  size_t pos = beg + len * chunkNum / numOfChunks;
  size_t seqNum = m.whichSequence(pos);
  size_t begDistance = pos - m.seqBeg(seqNum);
  size_t endDistance = m.padEnd(seqNum) - pos;
  return (begDistance < endDistance) ? seqNum : seqNum + 1;
}

}  // end namespace cbrc
#endif  // MULTISEQUENCE_HH
