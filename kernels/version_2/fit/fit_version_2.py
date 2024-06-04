from warnings import filterwarnings
filterwarnings('ignore')

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.model_selection import RepeatedKFold, train_test_split, GridSearchCV, KFold
from sklearn.compose import TransformedTargetRegressor
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error
from sklearn import preprocessing
from sklearn.gaussian_process.kernels import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import CSS4_COLORS as css

from sys import argv
import numpy as np
import pickle
import util
from util import bcolors as tcol
from tqdm import tqdm

plt.rcParams.update(util.my_rcParams)


    
if __name__=='__main__':
    """
        Plan of action:
            * split the kappa and Chi-squared set into test and training sets
            * construct a set of models 
            * use RepeatedKFold cross validation to find the best model
            * plot the result and save model to fole
    """
    rate_set = int(argv[1])
    rate_set_name = util.rate_set_names[rate_set]
    print(tcol.BOLD+tcol.OKGREEN+f"** Fitting For Rate Set: {rate_set} ** "+tcol.ENDC)
    
    ## Data import, feature engineering:
    data = util.read_data(rate_set)
    X = [[np.sqrt(item[0]), item[1]] for item in zip(data['1/kappa_r'],data['1/kappa_e'])]
    Y = data['chi_sq/ndf'].tolist()

    opt_starts = 5
    nrepeat, nsplit = 50, 5

    lscale_2d = (0.5, 0.5)
    lscale_1d =  0.2

    models = { 0 :GaussianProcessRegressor(kernel=RBF(length_scale=lscale_1d)            + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 1 :GaussianProcessRegressor(kernel=RBF(length_scale=lscale_1d)            + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 2 :GaussianProcessRegressor(kernel=Matern(nu=0.5, length_scale=lscale_1d) + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 3 :GaussianProcessRegressor(kernel=Matern(nu=1.5, length_scale=lscale_1d) + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 4 :GaussianProcessRegressor(kernel=Matern(nu=2.5, length_scale=lscale_1d) + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 5 :GaussianProcessRegressor(kernel=Matern(nu=0.5, length_scale=lscale_2d) + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 6 :GaussianProcessRegressor(kernel=Matern(nu=1.5, length_scale=lscale_2d) + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 7 :GaussianProcessRegressor(kernel=Matern(nu=2.5, length_scale=lscale_2d) + WhiteKernel(), n_restarts_optimizer=opt_starts)
             , 8 :GaussianProcessRegressor(kernel=Matern(nu=0.5, length_scale=lscale_1d) + RBF(length_scale=lscale_2d), n_restarts_optimizer=opt_starts)
             , 9 :GaussianProcessRegressor(kernel=Matern(nu=1.5, length_scale=lscale_1d) + RBF(length_scale=lscale_2d), n_restarts_optimizer=opt_starts)
             , 10:GaussianProcessRegressor(kernel=Matern(nu=2.5, length_scale=lscale_1d) + RBF(length_scale=lscale_2d), n_restarts_optimizer=opt_starts)
             , 11:GaussianProcessRegressor(kernel=RationalQuadratic(length_scale=lscale_1d,alpha=1) + WhiteKernel(), n_restarts_optimizer=opt_starts)}
    
    print(tcol.WARNING+tcol.BOLD+"* Model selection step: "+tcol.ENDC)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.10, train_size=0.90)
    scores, params = {}, {}
    gpr_params = {}
    best_model = None
    cv_min_score = 1e20
    
    for ik in models:
        scores[ik] = []
        params[ik] = []
        gpr_params[ik] = []
        model = models[ik]
        cross_validator = RepeatedKFold(n_splits=nsplit, n_repeats=nrepeat)
        pbar = tqdm(total=nrepeat*nsplit, colour='red')
        for ifold, (train_indices, test_indices) in enumerate(cross_validator.split(X_train)):
            pbar.update(1)
            x_train_set = [X_train[i] for i in train_indices]
            y_train_set = [Y_train[i] for i in train_indices]
            x_test_set  = [X_train[i] for i in test_indices]
            y_test_set  = [Y_train[i] for i in test_indices]
            model.fit(x_train_set, y_train_set)
            yhat = model.predict(x_test_set)
            score = mean_squared_error(y_true=y_test_set, y_pred=yhat, squared=False)
            if score < cv_min_score:
                cv_min_score = score
                best_model = model
        pbar.close()
        print(tcol.WARNING+f"Done with {ik} model family."+tcol.ENDC)
    
    yhat = model.predict(X_test)
    test_set_score = mean_squared_error(y_pred=yhat, y_true=Y_test, squared=False)
    print(tcol.HEADER+tcol.OKBLUE+f"Chosen model for rate set {rate_set} is:"+tcol.ENDC)
    print(tcol.BOLD+tcol.FAIL,best_model)
    print(f"with test set score of {test_set_score}", tcol.ENDC)

    ## Save information about the model in a log file and pickle the regressor
    with open(f"../../fit_results/chosen_mode_{rate_set}.pickle", 'wb') as f:
        pickle.dump(best_model, f)
    ## Plot:
    log_file = open(f"../../fit_results/Log_rset_{rate_set}.log",'w')
    log_file.write(f"Model number chosen with negative mean-squared-error: \n\t{cv_min_score:0.5e}\n")
    log_file.write(f"Model test set score: \n\t{test_set_score:0.5e}\n")
    log_file.write(f"Best Parameters: \n\t{best_model.kernel_.get_params()}\n")
    ## Test plot the model here:
    kapr = np.linspace(1, 15, 200, endpoint=True)
    kape = np.linspace(1, 15, 200, endpoint=True)
    grid_x,grid_y = np.meshgrid(1./kapr, 1./kape)
    plot_x,plot_y = np.meshgrid(kapr, kape)

    GRID_X = np.stack([grid_x.ravel(), grid_y.ravel()]).T
    PLOT_X = np.stack([plot_x.ravel(), plot_y.ravel()]).T
    
    
    fig1, ax1 = plt.subplots(1,1,figsize=(9,9))

    z              = best_model.predict(GRID_X)
    image          = np.reshape(z, grid_x.shape)
    minimum_chi_sq = min(z)
    print(minimum_chi_sq)
    indx_opt       = np.argmin(image)
    z_min, z_max   = 0, minimum_chi_sq+2.
    levels         = np.linspace(z_min, z_max, 5)
    result = PLOT_X[indx_opt]
    print_minimum_chi_sq = minimum_chi_sq
    log_file.write(f"{print_minimum_chi_sq:0.5e}\n")
    log_file.write(f"Minimum Chi-sq after the GPR: {print_minimum_chi_sq:0.5f}\n results: {result}")
    divider = make_axes_locatable(ax1)
    cax     = divider.append_axes('right', size='5%', pad=0.05)
    im      = ax1.pcolormesh(plot_x, plot_y, image, cmap='cividis_r', vmin=z_min, vmax=z_max, rasterized=True)
    cbar = fig1.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(r'$\chi^2/$n.d.f')
    cont   = ax1.contour(plot_x, plot_y, image, cmap='Reds_r', vmin=z_min, vmax=z_max, levels=levels)
    ax1.clabel(cont, inline=True, fontsize=15)
    ax1.set_ylabel('$\kappa_e$', fontsize=30)
    ax1.set_xlabel('$\kappa_r$', fontsize=30)
    ax1.set_title(f'{rate_set_name}')
    plt.show()
    fig1.savefig(f"../../fit_results/Test_Plot_Rate_Set_{rate_set}.png",dpi=150)
    log_file.close()
