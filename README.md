# olfactionInContextGLM

1 - Create data file for GLM:

```

% Indicar fn, que es el archivo ExperimentData, por ejemplo:
fn='/media/Data/SmartboxRecording_L3_20210910-184935_ExperimentData.mat';

% Elegir fnRaw, que es el nombre del archivo donde se exporta los datos en formato GLM, por ejemplo:
fnRaw='/media/Data/glmFits/glmDataL3_20210910-184935.mat';

% Correr:
rawDataForGLM(fn,fnRaw)

```



2 - Fit GLM to data:

```

% Ajustar todos los modelos a la lista neuronsToFit de neuronas a analizar

neuronsToFit=1:10; % e.g, fit 10 first neurons
fitModel='all';
fitGLM_lab(fnRaw,neuronsToFit,fitModel)
% FOR THIS fnRaw, RESULTS WILL BE SAVED IN /media/Data/glmFits/glmDataL3_20210910-184935_modelFit/

```



3 - Estimate the contribution of each variable to the fitted model

```

% Indicar fnBase, el nombre base de los archivos ajustados, en este ejemplo:
fnBase = '/media/Data/glmFits/glmDataL3_20210910-184935_modelFit/fit_Neuron';

% Calculamos las contribuciones
neuronsToFit=8; % e.g, calculate contribution of neuron 8
parameterContribution(fnBase,neuronsToFit)

```



4 - Remove variables with negative contributions and refit model:

```

% Correr:
checkAndCorrectModels(fnBase)

```



5 - Plot kernels obtained:

```

% Indicar neurona a graficar:
fnNeuronFit='/media/Data/Dropbox-labo/Dropbox/seba/Datos/Datos in vivo/glmFits/glmDataL3_20210910-184935_modelFit/fit_Neuron01_XOLIRSM.mat';
includeHP=1; % PARA INCLUIR HISTORY AND POP
plotFlag=[1 0]; % [1 0] significa agrupar variables (ej., olor 1 y olor 2) en un mismo grafico , [1 1] significa no agrupar
wKern=getKernels(fnNeuronFit,includeHP,plotFlag);

```
