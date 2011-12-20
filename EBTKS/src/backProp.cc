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
$RCSfile: backProp.cc,v $
$Revision: 1.4 $
$Author: jason $
$Date: 2004-01-19 15:38:15 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "backProp.h"
#include <assert.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MString.h"

using namespace std;


const int   BP_ANN::_LUT_LENGTH    = 1001;  
// Should be odd, as it will be used for a range of values centered around 0
const long  BP_ANN::_SEED          = 7366498;
const int   BP_ANN::_N_CYCLES      = 500;
const float BP_ANN::_LEARNING_RATE = 0.8;
const float BP_ANN::_MOMENTUM      = 0.3;
const float BP_ANN::_TEMPERATURE   = 1.0;
const float BP_ANN::_MAX_ERROR     = 0.0;
const float BP_ANN::_MAX_D_ERROR   = 1e-6;
const unsigned BP_ANN::_SHUFFLE_INTERVAL = 2*unsigned(MAXINT);
//const unsigned BP_ANN::_SHUFFLE_INTERVAL = 25;

BP_ANN::BP_ANN(const UnsignedArray& topology, Boolean verbose)
{
  _nLayers = _nInputNodes = _nOutputNodes = 0;
  _nNodesInLayer = NULL;
  _nWeightsForLayer = NULL;
  _node = NULL;
  _weight = NULL;
  _verbose = verbose;

  // Create network
  _create(topology);

  // Randomize network
  randomize(_SEED);

  // Set default parameters
  setDefaults();

  if (_verbose)
    save(cout, FALSE);
}

BP_ANN::BP_ANN(istream& IS, Boolean verbose)
{
  _nLayers = _nInputNodes = _nOutputNodes = 0;
  _nNodesInLayer = NULL;
  _nWeightsForLayer = NULL;
  _node = NULL;
  _weight = NULL;
  _verbose = verbose;

  load(IS);

  if (_verbose)
    save(cout, FALSE);
}

BP_ANN::~BP_ANN()
{
  _destroy();
}

Boolean 
BP_ANN::nInputNodes(unsigned n)
{
  if (n == _nInputNodes)
    return SUCCESS;

  if (!_nLayers) {
    cerr << "#Layers: " << _nLayers << endl;
    return FAILURE;
  }
  
  UnsignedArray topology(_nNodesInLayer, _nLayers);
  
  topology[(unsigned int)0] = n;

  _create(topology);
  randomize(_SEED);

  if (_verbose)
    save(cout, FALSE);

  return SUCCESS;
}

Boolean 
BP_ANN::nOutputNodes(unsigned n)
{
  if (n == _nOutputNodes)
    return SUCCESS;

  if (!_nLayers) {
    cerr << "#Layers: " << _nLayers << endl;
    return FAILURE;
  }
  
  UnsignedArray topology(_nNodesInLayer, _nLayers);
  
  topology[_nLayers - 1] = n;

  _create(topology);
  randomize(_SEED);

  if (_verbose)
    save(cout, FALSE);

  return SUCCESS;
}

Boolean 
BP_ANN::topology(const UnsignedArray& topology)
{
  _create(topology);
  randomize(_SEED);

  if (_verbose)
    save(cout, FALSE);

  return SUCCESS;
}

void
BP_ANN::randomize(long seed)
{
  srand48(seed ? seed : time(0));

  unsigned     layerCtr, i;
  unsigned     nNodes;
  unsigned     nWeights;
  WEIGHT *weightPtr;
  BPNODE *nodePtr;

  for (layerCtr = 0; layerCtr < _nLayers; layerCtr++) {

    nodePtr = _node[layerCtr];
    nNodes  = _nNodesInLayer[layerCtr];

    for (i = 0; i < nNodes; i++) {
      nodePtr->out = 0.0;
      nodePtr->bias = drand48() - 0.5;
      nodePtr->dBias = 0.0;
      nodePtr->delta = 0.0;
      nodePtr++;
    }
  }

  for (layerCtr = 1; layerCtr < _nLayers; layerCtr++) {

    weightPtr  = _weight[layerCtr];
    nWeights   = _nWeightsForLayer[layerCtr];

    for (i = 0; i < nWeights; i++) {
      weightPtr->value  = drand48() - 0.5;
      weightPtr->delta = 0.0;
      weightPtr++;
    }
  }
}

void
BP_ANN::setDefaults()
{
  _stopRequest  = FALSE;
  _learningRate = _LEARNING_RATE;
  _momentum     = _MOMENTUM;
  _nCycles      = _N_CYCLES;
  _cycleCtr     = 0;
  _nSamples     = 0;
  _maxError     = _MAX_ERROR;
  _maxDerror    = _MAX_D_ERROR;
  _shuffleInterval      = _SHUFFLE_INTERVAL;
  _sigmoidAtOutputLayer = TRUE;
  _softmaxAtOutputLayer = FALSE;

  _createLut(_TEMPERATURE);
}

void
BP_ANN::set(double alpha, double beta, double T, unsigned nCycles,
	    double maxError, double maxDerror, unsigned shuffleInterval, Boolean sigmoid,
	    Boolean softmax)
{
  _learningRate = alpha;
  _momentum     = beta;
  _nCycles      = nCycles;
  _maxError     = maxError;
  _maxDerror    = maxDerror;
  _shuffleInterval      = shuffleInterval;
  _sigmoidAtOutputLayer = sigmoid;
  _softmaxAtOutputLayer = softmax;

  _createLut(T);
}
/*
int
BP_ANN::train(TrainingSet& trainingSet, ErrorMonitor errorMonitor)
{
  if (initTraining(trainingSet.size()) != SUCCESS) {
    cerr << "Couldn't initialize training" << endl;
    return(FAILURE);
  }

  unsigned nOutputNodes = _nNodesInLayer[_nLayers - 1];
  unsigned totalOutputs = _nSamples*nOutputNodes;
 
  double previousMSE = MAXDOUBLE;
  unsigned shuffleInterval = _shuffleInterval;

  if (shuffleInterval) {
    cout << "Shuffling training set (" << trainingSet.size() << ")..." << flush;
    trainingSet.shuffle();
    cout << "Done" << endl;
  }

  TrainingSetIterator examplesIt(trainingSet);
  Example *example;

  for (unsigned cycleCtr = 0; !_stopRequest; cycleCtr++) {
    examplesIt.first();
    unsigned exampleCtr = 0;
    double   MSE = 0;
    while (example = examplesIt++) {
      train(example->input, example->target);
      MSE += _outputError.sum2();
      exampleCtr++;
    }  
    MSE /= totalOutputs;
    double dError = previousMSE - MSE;

    _stopRequest = 
      ((cycleCtr >= _nCycles) || 
       ((dError >= 0) && (dError <= _maxDerror)) ||
       (MSE <= _maxError));

    previousMSE = MSE;

    if (errorMonitor)
      (errorMonitor)(cycleCtr, MSE);
    
    if (shuffleInterval && !--shuffleInterval) {
      cout << "Shuffling training set..." << flush;
      trainingSet.shuffle();
      cout << "Done" << endl;
      shuffleInterval = _shuffleInterval;
      previousMSE = MAXDOUBLE;
    }
  }

  return(SUCCESS);
}
*/

int
BP_ANN::train(TrainingSet& trainingSet, ErrorMonitor errorMonitor)
{
  if (initTraining(trainingSet.size()) != SUCCESS) {
    cerr << "Couldn't initialize training" << endl;
    return(FAILURE);
  }

  int nextSample = 0;
  const Example *sample = &trainingSet[nextSample];
  while ((nextSample = train(sample->input, sample->target, errorMonitor)) >= 0)
    sample = &trainingSet[nextSample];
  
  return(SUCCESS);
}

int
BP_ANN::initTraining(unsigned nSamples)
{
  if (!_lut.size()) {
    cerr << "Error: LUT not created!" << endl;
    return(FAILURE);
  }

  if (_nLayers <= 0) {
    cerr << "Error: Invalid # layers (" << _nLayers << ")" << endl;
    return(FAILURE);
  }

  if (!nSamples) {
    cerr << "Error: #samples: " << nSamples << endl;
    return(FAILURE);
  }

  _cycleCtr    = 0;
  _nSamples    = nSamples; 
  _sampleCtr   = 0;
  _stopRequest = FALSE;

  return SUCCESS;
}

int
BP_ANN::train(const double *input, const double *target, ErrorMonitor errorMonitor)
{
  if (!_nSamples) {
    cerr << "Error: #samples: " << _nSamples << endl;
    return -1;
  }

  static double                previousMSE, MSE;
  static SimpleArray<unsigned> sampleList;
  static unsigned              shuffleInterval;
  static unsigned              totalOutputs;

  if (!_sampleCtr) {
    if (!_cycleCtr) { // Starting first cycle
      totalOutputs    = _nSamples*_nOutputNodes;
      shuffleInterval = 0; // Force shuffle at first cycle
      sampleList      = SimpleArray<unsigned>(0, 1, _nSamples - 1);
      MSE             = 0;
    }
    else {
      MSE /= totalOutputs;
      double dError = previousMSE - MSE;

      if (errorMonitor)
	(errorMonitor)(_cycleCtr, MSE);

      if (_stopRequest || 
	  (_cycleCtr >= _nCycles) || 
	  ((dError >= 0) && (dError <= _maxDerror)) ||
	  (MSE <= _maxError))
	return -1;

      previousMSE = MSE;
      MSE = 0;
    }
    
    if (!shuffleInterval--) {
      if (_verbose)
	cout << "Shuffling training set..." << flush;
      sampleList.shuffle();
      if (_verbose)
	cout << "Done" << endl;
      shuffleInterval = _shuffleInterval - 1;
      previousMSE = MAXDOUBLE;
    }

    _sampleCtr = _nSamples;
    _cycleCtr++;
  }
  
  _forward(input);
  _calculateDeltas(target);
  _adjustWeights();

  MSE += _outputError.sum2();

  /*
  cout << sampleList << "; " 
       << _sampleCtr << ", " 
       << _cycleCtr << ", " 
       << shuffleInterval << endl;
       */

  return int(sampleList[--_sampleCtr]);
}

void
BP_ANN::evaluate(const double *input, double *output)
{
  BPNODE *nodePtr = _node[_nLayers - 1];

  _forward(input);

  for (register unsigned i = _nOutputNodes; i; i--)
    *output++ = (nodePtr++)->out;
}

unsigned
BP_ANN::classify(const double *input, double *output)
{
  BPNODE *nodePtr = _node[_nLayers - 1];

  _forward(input);

  double   maxVal   = -MAXDOUBLE;
  unsigned node     = 0;

  for (register unsigned i = 0; i < _nOutputNodes; i++, nodePtr++) {
    if (nodePtr->out > maxVal) {
      maxVal = nodePtr->out;
      node   = i;
    }
    if (output)
      *output++ = nodePtr->out;
  }

  return node;
}

ostream&
BP_ANN::printWeights(ostream& theStream)
{
  for (unsigned i = 1; i < _nLayers; i++) {
    for (unsigned j = 0; j < _nWeightsForLayer[i]; j++)
      theStream << _weight[i][j].value << " ";
    theStream << endl;
  }

  return theStream;
}

ostream&
BP_ANN::printNodes(ostream& theStream)
{
  for (unsigned i = 0; i < _nLayers; i++) {
    for (unsigned j = 0; j < _nNodesInLayer[i]; j++)
      theStream << _node[i][j].out << " ";
    theStream << endl;
  }

  return theStream;
}


int
BP_ANN::save(ostream& OS, Boolean includeContents)
{
  if (!OS)
    return FAILURE;

  unsigned layerCtr;
  OS << "learning_rate:     " << _learningRate << endl
     << "momentum:          " << _momentum << endl
     << "temperature:       " << _temperature << endl
     << "num_of_cycles:     " << _nCycles << endl
     << "max_error:         " << _maxError << endl
     << "max_d_error:       " << _maxDerror << endl
     << "shuffle_interval:  " << _shuffleInterval << endl
     << "layers:            " << _nLayers;
  for (layerCtr = 0; layerCtr < _nLayers; layerCtr++)
    OS << " " << _nNodesInLayer[layerCtr];
  OS << endl;

  if (includeContents) {
    OS << "contents:" << endl;
    for (layerCtr = 1; layerCtr < _nLayers; layerCtr++) {
      for (unsigned i = 0; i < _nNodesInLayer[layerCtr]; i++)
	OS << _node[layerCtr][i].bias << " "
	   << _node[layerCtr][i].dBias << " "
	   << _node[layerCtr][i].delta << " "
	   << _node[layerCtr][i].out << endl;
      OS << endl;
    }
    
    for (layerCtr = 1; layerCtr < _nLayers; layerCtr++) {
      for (unsigned i = 0; i < _nWeightsForLayer[layerCtr]; i++)
	OS << _weight[layerCtr][i].value << " "
	   << _weight[layerCtr][i].delta << endl;
      OS << endl;
    }
  }

  return SUCCESS;
}

int
BP_ANN::load(istream& IS)
{
  // Initialize with default values
  setDefaults();
  double        temperature = _TEMPERATURE;
  long          seed        = _SEED;
  UnsignedArray tplgy;
  Boolean       contentsRead = FALSE;

  MString key;
  while (IS >> key) {
    if (key.contains("randomize"))
      IS >> seed;
    else if (key.contains("learning_rate"))
      IS >> _learningRate;
    else if (key.contains("momentum"))
      IS >> _momentum;
    else if (key.contains("temperature"))
      IS >> temperature;
    else if (key.contains("num_of_cycles"))
      IS >> _nCycles;
    else if (key.contains("max_error"))
      IS >> _maxError;
    else if (key.contains("max_d_error"))
      IS >> _maxDerror;
    else if (key.contains("shuffle_interval"))
      IS >> _shuffleInterval;
    else if (key.contains("layers")) {
      unsigned nLayers;
      IS >> nLayers;
      tplgy.newSize(nLayers);
      for (unsigned layerCtr = 0; layerCtr < nLayers; layerCtr++)
	IS >> tplgy[layerCtr];
      _create(tplgy);
    }
    else if (key.contains("contents")) {
      unsigned i, layerCtr;
      // Read node values
      for (layerCtr = 1; layerCtr < _nLayers; layerCtr++)
        for (i = 0; i < _nNodesInLayer[layerCtr]; i++)
	  IS >> _node[layerCtr][i].bias
	     >> _node[layerCtr][i].dBias
	     >> _node[layerCtr][i].delta
	     >> _node[layerCtr][i].out;
      
      // Read weights
      for (layerCtr = 1; layerCtr < _nLayers; layerCtr++)
        for (i = 0; i < _nWeightsForLayer[layerCtr]; i++)
	  IS >> _weight[layerCtr][i].value
	     >> _weight[layerCtr][i].delta;

      contentsRead = TRUE;
    }
  }

  _createLut(temperature);

  if (!contentsRead)
    randomize(seed);

  return(SUCCESS);
}

void
BP_ANN::_forward(const double *input)
{
  BPNODE *nodePtr, *prevLayerNodePtr;
  BPNODE *layerPtr, *prevLayerPtr;
  WEIGHT *weightPtr;
  double  sum;
  int          layerCtr;
  unsigned     nNodes, nNodesInPrevLayer, nodeCtr;
  unsigned     i;
  
  nNodes     = _nNodesInLayer[0];
  nodePtr    = _node[0];
  const double *inputPtr   = input;
  for (i = 0; i < nNodes; i++)
    (nodePtr++)->out = *inputPtr++;

  prevLayerPtr      = _node[0];
  nNodesInPrevLayer = _nNodesInLayer[0];

  double  *lut          = _lut.contents();
  int      lastLutIndex = _lut.size() - 1;

  // Hidden layers
  int lastHiddenLayer = _nLayers - 1;
  for (layerCtr = 1; layerCtr < lastHiddenLayer; layerCtr++) {

    weightPtr    = _weight[layerCtr];
    nodePtr      = layerPtr = _node[layerCtr];
    nNodes       = _nNodesInLayer[layerCtr];

    for (nodeCtr = 0; nodeCtr < nNodes; nodeCtr++) {

      prevLayerNodePtr = prevLayerPtr;
      sum = 0.0;

      for (i = 0; i < nNodesInPrevLayer; i++)
	sum += (weightPtr++)->value * (prevLayerNodePtr++)->out;

      int lutIndex = (int) rint((sum + nodePtr->bias + 5)/_lutStep);
      if (lutIndex < 0)
	lutIndex = 0;
      else if (lutIndex > lastLutIndex)
	lutIndex = lastLutIndex;

      (nodePtr++)->out = lut[lutIndex];
    }

    prevLayerPtr = layerPtr;
    nNodesInPrevLayer = nNodes;
  }

  // Output layer

  layerCtr  = _nLayers - 1;
  weightPtr = _weight[layerCtr];
  nodePtr   = _node[layerCtr];
  nNodes    = _nNodesInLayer[layerCtr];
  double outputSum = 0.0;
  
  for (nodeCtr = 0; nodeCtr < nNodes; nodeCtr++, nodePtr++) {
    prevLayerNodePtr = prevLayerPtr;
    sum = 0.0;

    // Disconnect one node if softmax is used
    if (!_softmaxAtOutputLayer || nodeCtr)
      for (i = 0; i < nNodesInPrevLayer; i++) {
	sum += (weightPtr++)->value * (prevLayerNodePtr++)->out;
      }
    
    if (_sigmoidAtOutputLayer) {
      int lutIndex = (int) rint((sum + nodePtr->bias + 5)/_lutStep);
      if (lutIndex < 0)
	lutIndex = 0;
      else if (lutIndex > lastLutIndex)
	lutIndex = lastLutIndex;
      
      nodePtr->out = lut[lutIndex];
    }
    else
      nodePtr->out = sum + nodePtr->bias;

    if (_softmaxAtOutputLayer) {
      outputSum += (nodePtr->out = exp(nodePtr->out));
      //      if (*input)
      //	cout << nodePtr->out << ":" << outputSum << " " << flush;
    }
  }
  
  if (_softmaxAtOutputLayer) {
    nodePtr = _node[layerCtr];
    for (nodeCtr = 0; nodeCtr < nNodes; nodeCtr++, nodePtr++)
      nodePtr->out /= outputSum;
  }
}

void
BP_ANN::_calculateDeltas(const double *target)
{
  BPNODE *nextLayerNodePtr;
  BPNODE *layerPtr, *nextLayerPtr;
  WEIGHT *weightLayerPtr;
  WEIGHT *weightPtr;
  double  nodeValue;
  double  sum;
  int     outputLayerIndex;
  int     layerCtr;
  int     nNodesInNextLayer;

/*
 * Calculate deltas for the output layer
 */

  outputLayerIndex = _nLayers - 1;

  double *outputErrorPtr  = _outputError.contents();
  BPNODE *nodePtr         = _node[outputLayerIndex];
  int nNodes              = _nNodesInLayer[outputLayerIndex];

  unsigned i;
  if (_sigmoidAtOutputLayer)
    for (i = nNodes; i != 0; i--) {
      nodeValue = nodePtr->out;
      *outputErrorPtr = *target++ - nodeValue;
      (nodePtr++)->delta = (*outputErrorPtr++)*nodeValue*(1.0 - nodeValue);
    }
  else
    for (i = nNodes; i != 0; i--)
      (nodePtr++)->delta = *outputErrorPtr++ = *target++ - nodePtr->out;
  
/*
 * Calculate deltas for the hidden layer(s)
 */

  nextLayerPtr = _node[outputLayerIndex];
  nNodesInNextLayer = _nNodesInLayer[outputLayerIndex];

  for (layerCtr = outputLayerIndex - 1; layerCtr > 0; layerCtr--) {
    nodePtr        = layerPtr = _node[layerCtr];
    weightLayerPtr = _weight[layerCtr + 1];
    nNodes         = _nNodesInLayer[layerCtr];

    for (register unsigned i = nNodes; i != 0; i--) {
      weightPtr        = weightLayerPtr++;
      nextLayerNodePtr = nextLayerPtr;

      sum = 0.0;
      for (register unsigned j = nNodesInNextLayer; j != 0; j--) {
	sum += weightPtr->value * nextLayerNodePtr->delta;

	nextLayerNodePtr++;
	weightPtr += nNodes;
      }
      nodePtr->delta = nodePtr->out * (1 - nodePtr->out) * sum;
      nodePtr++;
    }

    nNodesInNextLayer = nNodes;
    nextLayerPtr = layerPtr;
  }    
}

void 
BP_ANN::_create(const UnsignedArray& topology)
{
  _destroy();

  _nLayers = topology.size();

  if (!_nLayers) {
    cerr << "# network layers is zero!" << endl;
    exit(EXIT_FAILURE);
  }

  if ((_nNodesInLayer = (unsigned *) malloc(_nLayers*sizeof(unsigned))) == NULL)
    assert(0);

  unsigned i;
  for (i = 0; i < _nLayers; i++)
    _nNodesInLayer[i] = topology[i];

  _nInputNodes  = topology[(unsigned int)0];
  _nOutputNodes = topology[_nLayers - 1];

  if ((_nWeightsForLayer = (unsigned *) malloc(_nLayers*sizeof(unsigned))) == NULL)
    assert(0);
  for (i = 1; i < _nLayers; i++)
    _nWeightsForLayer[i] = _nNodesInLayer[i - 1]*_nNodesInLayer[i];

/*
 * Allocate node outputs for all layers
 */

  if ((_node = (BPNODE **) malloc(_nLayers*sizeof(BPNODE *))) == NULL)
    assert(0);
  for (i = 0; i < _nLayers; i++)
    if ((_node[i] = (BPNODE *) calloc(_nNodesInLayer[i], sizeof(BPNODE))) == NULL)
      assert(0);

/*
 * Allocate weights for all layers except the input layer
 */

  if ((_weight = (WEIGHT **) malloc(_nLayers*sizeof(WEIGHT *))) == NULL)
    assert(0);
  for (i = 1; i < _nLayers; i++)
    if ((_weight[i] = 
	 (WEIGHT *) calloc(_nWeightsForLayer[i], sizeof(WEIGHT))) == NULL)
      assert(0);

/*
 * Allocate memory for the outputErrors
 */

  _outputError.newSize(_nOutputNodes);
}

void
BP_ANN::_destroy()
{
  if (_nNodesInLayer != NULL) {
    free((unsigned *) _nNodesInLayer);
    _nNodesInLayer = NULL;
  }
  if (_nWeightsForLayer != NULL) {
    free((unsigned *) _nWeightsForLayer);
    _nWeightsForLayer = NULL;
  }
  if (_node != NULL) {
    for (unsigned i = 0; i < _nLayers; i++)
      free((BPNODE *) _node[i]);
    free((BPNODE **) _node);
    _node = NULL;
  }

  if (_weight != NULL) {
    for (unsigned i = 1; i < _nLayers; i++)
      free((WEIGHT *) _weight[i]);
    free((WEIGHT **) _weight);
    _weight = NULL;
  }

  _nLayers = 0;
}

void
BP_ANN::_createLut(double T)
{
  if ((_lut.size() == _LUT_LENGTH) && (_temperature == T))
    return;

  _lut.newSize(_LUT_LENGTH);

  _lutStep = 10.0/(_LUT_LENGTH - 1);
  for (register unsigned i = 0; i < _LUT_LENGTH; i++)
    _lut[i] = 1.0/(1.0 + exp(-(_lutStep*i - 5.0)/T));

  _temperature = T;
}
    
void
BP_ANN::_adjustWeights()
{
  BPNODE *layerPtr, *prevLayerPtr;
  BPNODE *nodePtr, *prevLayerNodePtr;
  WEIGHT *weightPtr;
  unsigned     layerCtr;
  unsigned     nNodes, nNodesInPrevLayer;
  unsigned     nodeCtr;
  unsigned     i;

  prevLayerPtr      = _node[0];
  nNodesInPrevLayer = _nNodesInLayer[0];

  for (layerCtr = 1; layerCtr < _nLayers; layerCtr++) {
    
    weightPtr  = _weight[layerCtr];
    nodePtr    = layerPtr = _node[layerCtr];
    nNodes     = _nNodesInLayer[layerCtr];

    for (nodeCtr = 0; nodeCtr < nNodes; nodeCtr++) {

      prevLayerNodePtr = prevLayerPtr;

      for (i = 0; i < nNodesInPrevLayer; i++) {
	weightPtr->delta = 
	  _learningRate * nodePtr->delta * prevLayerNodePtr->out +
	    _momentum * weightPtr->delta;
	weightPtr->value += weightPtr->delta;

	weightPtr++;
	prevLayerNodePtr++;
      }
      
      nodePtr->dBias = _learningRate*nodePtr->delta + _momentum*nodePtr->dBias;
      nodePtr->bias += nodePtr->dBias;
      
      nodePtr++;
    }

    prevLayerPtr      = layerPtr;
    nNodesInPrevLayer = nNodes;
  }
}
