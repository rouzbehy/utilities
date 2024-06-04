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
    print(util.bcolors.BOLD+util.bcolors.OKGREEN+f"Fitting For Rate Set: {rate_set}"+util.bcolors.ENDC)
    
    ## Data import, feature engineering:
    data = util.read_data(rate_set)
    X = [[item[0], item[1]] for item in zip(data['1/kappa_r'],data['1/kappa_e'])]
    Y = data['chi_sq/ndf'].tolist()
    ## Define the parameter grid to do the grid search, optimize each kernel
    lbounds = (1e-5,1e5)
    noise_bounds = (1e-4,1e0)
    

    ## Prep for model selection
    print(util.bcolors.WARNING+util.bcolors.BOLD+"\t Model selection step: "+util.bcolors.ENDC)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.30, train_size=0.70, random_state=0)
    ## the fitting:
    scaler = preprocessing.StandardScaler()
    scaler.fit(X_train)
    scaler.transform(X_train)
    cross_validator = RepeatedKFold(n_repeats=5, n_splits=4)
    for train_indices, test_indices in cross_validator.split(X_train):
        #print(f"\t train: {train_indices}\n\t test: {test_indices}")
        #print(train_indices, test_indices)
        kernel = RBF(length_scale=np.ones_like(X_train[0]), length_scale_bounds=lbounds) + WhiteKernel(noise_level_bounds=noise_bounds)
        model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=2)
        x_train_set = [X_train[i] for i in train_indices]
        y_train_set = [Y_train[i] for i in train_indices]
        x_test_set  = [X_train[i] for i in test_indices]
        y_test_set  = [Y_train[i] for i in test_indices]
        print(model)
        model.fit(x_train_set, y_train_set)
        print(model.kernel_)
        yhat = model.predict(x_test_set)
        err = mean_squared_error(y_test_set, yhat)
        print(err)
    # print(util.bcolors.WARNING+util.bcolors.BOLD, search.best_params_, search.best_score_, util.bcolors.ENDC)
    # model = search.best_estimator_
    # test_score = model.score(X_test, Y_test)
    # print(util.bcolors.FAIL+util.bcolors.BOLD, "test score: ", test_score, util.bcolors.ENDC)

    # ## Save information about the model in a log file and pickle the regressor
    # with open(f"../../fit_results/chosen_mode_{rate_set}.pickle", 'wb') as f:
    #     pickle.dump(model, f)

    # log_file = open(f"../../fit_results/Log_rset_{rate_set}.log",'w')
    # log_file.write(f"Model number chosen with negative mean-squared-error: \n\t{search.best_score_:0.5e}\n")
    # log_file.write(f"Model test set score: \n\t{test_score:0.5e}\n")
    # log_file.write(f"Best Parameters: \n\t{search.best_params_}\n")
    # log_file.write(f"Best Estimator: \n\t{search.best_estimator_.get_params()}\n")
    

    # ## Test plot the model here:
    # kapr = np.linspace(1, 15, 400, endpoint=True)
    # kape = np.linspace(1, 15, 400, endpoint=True)
    # grid_x,grid_y = np.meshgrid(1/kapr, 1/kape)
    # plot_x,plot_y = np.meshgrid(kapr, kape)

    # GRID_X = np.stack([grid_x.ravel(), grid_y.ravel()]).T
    # PLOT_X = np.stack([plot_x.ravel(), plot_y.ravel()]).T
    # #fig1, ax1 = plt.subplots(1,1,figsize=(7.5,7.5))

    # z              = model.predict(GRID_X)
    # image          = np.reshape(z, grid_x.shape)
    # minimum_chi_sq = min(z)
    # indx_opt       = np.argmin(image)
    # z_min, z_max   = 0, minimum_chi_sq+2.
    # levels         = np.linspace(z_min, z_max, 5)
    # result = PLOT_X[indx_opt]
    # print_minimum_chi_sq = minimum_chi_sq
    # log_file.write(f"{print_minimum_chi_sq:0.5e}\n")
    # log_file.write(f"Minimum Chi-sq after the GPR: {print_minimum_chi_sq:0.5f}\n results: {result}")
    # #divider = make_axes_locatable(ax1)
    # #cax     = divider.append_axes('right', size='5%', pad=0.05)
    # #im      = ax1.pcolormesh(plot_x, plot_y, image, cmap='cividis_r', vmin=z_min, vmax=z_max, rasterized=True)
    # #cbar = fig1.colorbar(im, cax=cax, orientation='vertical')
    # #cbar.set_label(r'$\chi^2/$n.d.f')
    # #cont   = ax1.contour(plot_x, plot_y, image, cmap='Reds_r', vmin=z_min, vmax=z_max, levels=levels)
    # #ax1.clabel(cont, inline=True, fontsize=15)
    # #ax1.set_ylabel('$\kappa_e$', fontsize=30)
    # #ax1.set_xlabel('$\kappa_r$', fontsize=30)
    # #ax1.set_title(f'{rate_set_name}')
    # #fig1.savefig(f"../../chisquared_work/Test_Plot_Rate_Set_{rate_set}_inverse_kappas_vs_chisq.png",dpi=150)
    # log_file.close()
    # #plt.show()
