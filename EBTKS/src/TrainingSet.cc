/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: TrainingSet.cc,v $
$Revision: 1.3 $
$Author: stever $
$Date: 2003-11-17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "TrainingSet.h"
#include <assert.h>
#include <stdlib.h>

using namespace std;



const int DEFAULT_TRAINING_SET_SIZE = 250;

/***************
 * Example class
 ***************/

//
// Constructors/destructor
//

Example::Example(unsigned nInputs, unsigned nTargets)
: input(0.0, nInputs),
  target(0.0, nTargets)
{}

Example::Example(unsigned lbl, const DblArray& inpt, const DblArray& trgt)
: input(inpt),
  target(trgt)
{
  label = lbl;
}

Example::~Example()
{}

ostream&
operator << (ostream& theStream, Example& example)
{
  theStream << example.label << ": " << example.input << " -> " << example.target;
  return theStream;
}

/*******************
 * TrainingSet class
 *******************/

TrainingSet::TrainingSet()
: OrderedCltn(DEFAULT_TRAINING_SET_SIZE)
{}

TrainingSet::TrainingSet(unsigned nExamples, unsigned nInputs, unsigned nTargets,
			 double minTarget, double maxTarget)
: OrderedCltn(nExamples)
{
  set(nInputs, nTargets, minTarget, maxTarget);
}

TrainingSet::~TrainingSet()
{
  removeAll();
}

//
// Set functions
//

void 
TrainingSet::set(unsigned nInputs, unsigned nTargets, double minTarget, double maxTarget)
{
  _nInputs   = nInputs;
  _nTargets  = nTargets;
  _minTarget = minTarget;
  _maxTarget = maxTarget;
}

//
// Adding examples
//

void 
TrainingSet::add(unsigned node, const double *inputValues)
{
  DblArray targetValues(_minTarget, _nTargets);
  targetValues[node] = _maxTarget;

  Example *newExample = new Example(node, DblArray(inputValues, _nInputs), targetValues);
  assert(newExample);

  //cout << "Adding example: " << *newExample << endl;

  OrderedCltn::add(newExample);
}

void 
TrainingSet::add(const double *iValues, const double *tValues)
{
  DblArray inputValues(iValues, _nInputs);
  DblArray targetValues(tValues, _nTargets);

  unsigned node;
  targetValues.max(&node);

  Example *newExample = new Example(node, inputValues, targetValues);
  assert(newExample);

  //cout << "Adding example: " << *newExample << endl;

  OrderedCltn::add(newExample);
}

void 
TrainingSet::add(unsigned class1, const double *inputValues1,
		 unsigned class2, const double *inputValues2, double pct)
{
  double targetRange   = _maxTarget - _minTarget;
  double pctComplement = 1.0 - pct;

  DblArray targetValues(_minTarget, _nTargets);
  targetValues[class1] = pct*targetRange + _minTarget;
  targetValues[class2] = pctComplement*targetRange + _minTarget;

  DblArray inputValues(_nInputs);
  for (unsigned i = 0; i < _nInputs; i++)
    inputValues[i] = pct*inputValues1[i] + pctComplement*inputValues2[i];

  Example *newExample = 
    new Example((pct < 0.5) ? class1 : class2, inputValues, targetValues);
  assert(newExample);

  //cout << "Adding example: " << *newExample << endl;

  OrderedCltn::add(newExample);
}

void
TrainingSet::add(const Array<SimpleArray<unsigned> >& mixtures, unsigned nMixtures)
{
  if (!nMixtures)
    return;

  // Setup mixture percentiles

  DblArray mixturePct(nMixtures);
  for (unsigned m = 0; m < nMixtures; m++)
    mixturePct[m] = double(m + 1)/(nMixtures + 1);

  // Create a separate TrainingSet for the new examples
  
  TrainingSet newExamples(uEndIndex, _nInputs, _nTargets, _minTarget, _maxTarget);

  // Cycle through all mixtures

  TrainingSetIterator class1it(*this);
  TrainingSetIterator class2it(*this);

  unsigned nMixedClasses = mixtures.size();
  for (unsigned class1 = 0; class1 < nMixedClasses; class1++) {
    const SimpleArray<unsigned>& combArray = mixtures[class1];
    unsigned nCombinations = combArray.size();
    for (unsigned comb = 0; comb < nCombinations; comb++) {
      unsigned class2 = combArray[comb];

      class1it.first();
      class2it.first();
      Example *class1ex, *class2ex;
      Boolean done = FALSE;

      while (!done) {
	// Search for next example of class 1

	while ((class1ex = class1it++) && (class1ex->label != class1));
	if (!class1ex)
	  done = TRUE;
	else {
	  // Search for next example of class 2

	  while ((class2ex = class2it++) && (class2ex->label != class2));
	  if (!class2ex)
	    done = TRUE;
	  else {
	    //cout << "Creating mixture for " << class1 << " and " << class2 << endl;

	    // Create new examples from the class1 and class2 examples

	    for (unsigned m = 0; m < nMixtures; m++)
	      newExamples.add(class1, class1ex->input, class2, class2ex->input, 
			      mixturePct[m]);
	  }
	}
      }
    }
  }

  this->addAllLast(newExamples);

  newExamples.uEndIndex = 0; // Prevent all added examples from being deleted
}

void
TrainingSet::add(const Array<SimpleArray<unsigned> >& mixtures, unsigned nMixtures, 
		 unsigned nCopies)
{
  if (!nMixtures)
    return;

  // Setup mixture percentiles

  DblArray mixturePct(nMixtures);
  for (unsigned m = 0; m < nMixtures; m++)
    mixturePct[m] = double(m + 1)/(nMixtures + 1);

  // Create a separate TrainingSet for the new examples
  
  TrainingSet newExamples(uEndIndex, _nInputs, _nTargets, _minTarget, _maxTarget);

  // Cycle through all mixtures

  TrainingSetIterator classIt(*this);

  unsigned nMixedClasses = mixtures.size();
  for (unsigned class1 = 0; class1 < nMixedClasses; class1++) {
    DblArray class1mean(0.0, _nInputs);

    unsigned exampleCtr = 0;
    classIt.first();
    Example *classEx;
    while (classEx = classIt++)
      if (classEx->label == class1) {
	class1mean += classEx->input;
	exampleCtr++;
      }

    if (exampleCtr > 0) {
      class1mean /= exampleCtr;
    
      const SimpleArray<unsigned>& combArray = mixtures[class1];
      unsigned nCombinations = combArray.size();
      for (unsigned comb = 0; comb < nCombinations; comb++) {
	unsigned class2 = combArray[comb];
	DblArray class2mean(0.0, _nInputs);

	unsigned exampleCtr = 0;
	classIt.first();
	Example *classEx;
	while (classEx = classIt++)
	  if (classEx->label == class2) {
	    class2mean += classEx->input;
	    exampleCtr++;
	  }

	if (exampleCtr > 0) {
	  class2mean /= exampleCtr;

	  for (unsigned n = nCopies; n != 0; n--)
	    for (unsigned m = 0; m < nMixtures; m++)
	      newExamples.add(class1, class1mean, class2, class2mean, mixturePct[m]);
	}
      }
    }
  }
  
  this->addAllLast(newExamples);

  newExamples.uEndIndex = 0; // Prevent all added examples from being deleted
}

//
// Removing examples
//

void
TrainingSet::removeAll()
{
  ocIterator exampleIt(*this);
  Example   *example;

  while (example = (Example *) exampleIt++)
    delete example;

  uEndIndex = 0;
}

//
// Misc functions
//

void
TrainingSet::shuffle()
{
  for (unsigned i = 0; i < uEndIndex; i++) {
    unsigned j = unsigned(drand48()*uEndIndex);
    if (i != j) {
      void *tempPtr = Contents[i];
      Contents[i] = Contents[j];
      Contents[j] = tempPtr;
    }
  }
}

ostream&
TrainingSet::print(ostream& OS) const
{
  TrainingSetIterator theIt(*this);
  Example             *example;
  while (example = theIt++)
    OS << example->label << ": " << example->input << " -> " << example->target << endl;

  return OS;
}

/****************************
 * TrainingSet iterator class
 ****************************/

Example *
TrainingSetIterator::operator ++() 
{ 
  if (++_exampleIndex >= _trainingSet->uEndIndex) 
    return(0); 

  return (Example *) _trainingSet->Contents[_exampleIndex];
}

Example *
TrainingSetIterator::operator ++(int) 
{ 
  if (_exampleIndex >= _trainingSet->uEndIndex) 
    return(0); 

  return (Example *) _trainingSet->Contents[_exampleIndex++];
}

Example *
TrainingSetIterator::operator --() 
{ 
  if ((int) --_exampleIndex < 0) 
    return(0);
  return (Example *) _trainingSet->Contents[_exampleIndex];
}

Example *
TrainingSetIterator::operator --(int) 
{ 
  if ((int) _exampleIndex < 0) 
    return(0);
  return (Example *) _trainingSet->Contents[_exampleIndex--];
}

