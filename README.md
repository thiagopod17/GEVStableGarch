#   GEVStableGarch

R package for ARMA-GARCH or ARMA-APARCH estimation with GEV (Generalized Extreme Value) or Stable distributions

#   Things to do
    - Write a better deffinition of stable distribution, see def. 1.6 Nolan - chapter 1.
    - Change name of variables that are functions or objects in R, such as: gamma, data, etc.
    gamma need to be changed in function GSgarch.Fit
    - Test the 0-parametrization indicated from Nolan to see if it performs better. Then,
    make sure the definition of the stable distribution is correct on the paper.
    - How to cite the code I used from package fGarch: Sim and garchSpec functions.
    - Include the gev restriction of stationarity into the model and test it.
    - Test GEV density to see why qsi < 0 is somewhat frequent in real data estimation.     Notice that the Zhao paper has qsi > 0.
    - Install the package "stable" from Nolan, J.P. and test the function GSgarch.Fit.
    - Test the estimation using the "sqp.restriction" algorithm.
    - Test the stable.moment.stationarity function using the particular cases provided 
    in Mittinik et al (2002).
    - Implementar overload do metodo print para objetos da classe fGEVSTABLEGARCH
    - Prediction of conditional variance for GARCH/APARCH models.
    -Separar funções ARMA com GEV e estavel das funções ARMA-GARCH.
    - Modify function names to make the package more user friendly.
    - Cant we define the ARMA-APARCH model with t-Student distribution with infinite 
    variance? degrees of freedom > 0 and not > 2.
    - Correct ARMA-stable estimation
      1.1 get.start para stable sem pegar variancia. pegar mean(abs(data-Mean))
      1.2 garch.Dist, fazer diferente. Calcular density stable com parametros z 
      (data) e hh(para sigma da estavel)
      1.3 definir valor de N para filter.Arma
      1.4 Ver sobre valores de inicio, superiores e inferiores do sigma da estavel, da media do processo e dos coeficientes arma.

#   Important notes about the modifications

    We will make the interface of this function more similar to the garchFit function
    from package fGarch to make it more user friendly. 
    The input parameters for this version (march/2015) are:

    - 25 Feb, 2015.
    Parameter input modification
    OK. algorithm: a string parameter that determines the algorithm used for maximum   likelihood estimation.
    OK. cond.dist: name of the conditional distribution, one of gev, stable, norm, std, sstd
    OK. control: control parameters, the same as used for the functions from nlminb, and 'bfgs' and 'Nelder-Mead' from optim.
    OK. data: The dataset to be estimated.
    formula: formula object describing the mean and variance equation of the ARMA-GARCH/APARCH model.
    OK. intercept: this flag determines if the parameter for the mean will be estimated or not
    OK. print.Result (Padrao eh TRUE): A boolean variable specifying whether or not the user wants to print the results after the function calling.
    OK. get.res: (NAO VAMOS TER MAIS ESSA VARIAVEL)
    OK. GSstable.tol e GStol: (CONFIGURAR NO INICIO DA FUNCAO, NAO MAIS NECESSARIA AQUI)
    APARCH: Vamos tirar pois usaremos da formula.
    Before Commiting to github: We saw that the estimated parameters from both GEVStableGarch package from CRAN and from our current version are the same. On the other hand, we saw that the estimated parameters from macbook differ slightly from the windows version. Our goal now is to investigate the filtering process inside the GSGarch.Fit function to make the estimated parameters more similar to the ones from package fGarch. 
    
    - 27 feb, 2015, right before commiting to Github
    Using garch11Fit function from Wurtz (2006) to estimate pure garch(1,1) model with conditional normal distribution This function estimate the garch(1,1)-include.mean-norm-dem2gbp
    
    - The results are exactly the same as in the Code Snippet 2 presented in the papper Wurtz et al. (2006) The function garch11Fit works better if start the conditional  variance with 'var(x)'. Mehoramos muito minha funcao quando para a estimacao do garch(1,1). Fiz isso retirando o filtro do aparch e recolocando o filtro do wuertz, que funciona para o garch11.
 
    - When I put the filter from garch11Fit function inside my GSgarch.Fit the results were exactly the same. Therefore, thats is our starting pointing. Now, I am almost done because my filter function for pure APARCH model is matching exactly the filter function from garch11Fit.   
    