# olfactionGLM

This is the repository for the MATLAB code used to fit GLM models in the paper LINK

## Acknowledgements 
The code is an adaptation of the [neuroGLM](https://github.com/pillowlab/neuroGLM) framework, for regressing trial-based spike train data using a Generalized Linear Model (GLM), introduced by [Park et al., 2014](https://pillowlab.princeton.edu/pubs/abs_ParkI_NN14.html).

It also makes use of tools from the [NCCLABCODE](https://github.com/pillowlab/DRD/tree/ddb2683d95fa4887156204ff472028ddd1dbb44b/ncclabcode) repository.

The cross-validation procedure for model selection is inspired by [Hardcastle et al., 2017](https://www.cell.com/neuron/fulltext/S0896-6273(17)30237-4) (implemented in [ln-model-of-mec-neurons](https://github.com/GiocomoLab/ln-model-of-mec-neurons))

## General
We provide a file piriformData.mat that contains data for 182 trials performed by an animal expert in the task. The data consists of the 10-ms time-binned spike count activity for 12 neurons, along with all the other covariates recorded during the task. These covariates (and their labels) are: Position in virtual corridor (X), Odorant stimulation (O), Licking response (L), Inhalation onset (I), Reward consumption (R), Pre-GO time window (G) and Contextual modulation of odor response (M).

## Reference
Paper citation LINK

## Instructions 

1 - Fit GLM to data:

```

% Let's use run the model selection and fitting procedure
fn='test dataset/piriformData.mat'; % Indicate file name with data
neuronsToFit=1:10; % e.g, fit 10 first neurons
fitModel='all'; % Indicate if a specific model is desired (eg., for a position model indicate 'X'), otherwise indicate 'all' to run thorugh the model selection procedure.
fitGLM_lab(fn,neuronsToFit,fitModel)
% FOR THIS fn, RESULTS FOR EACH NEURON WILL BE SAVED IN THE test dataset/piriformData_modelFit/ FOLDER

```

2 - Estimate the contribution of each variable to the fitted model

```
fnBase = 'test dataset/piriformData_modelFit/fit_Neuron'; % Indicate basename for neuron fit files
neuronsToFit=8; % e.g, calculate contribution of neuron 8
parameterContribution(fnBase,neuronsToFit)

```

3 - Plot kernels obtained:

```

% Indicate neuron fit file to plot:
fnNeuronFit='test dataset/piriformData_modelFit/fit_Neuron/fit_Neuron01_OLISM.mat'; % For neuron 1 of this dataset, the model selected is the model containing OLISM variables
includeHP=0; % Do not include history or population-coupling kernels
plotFlag=[1 0]; % 
wKern=getKernels(fnNeuronFit,includeHP,plotFlag);

```
