#ARIMA Model
import warnings
from pandas import read_csv
from pandas import datetime
from statsmodels.tsa.arima_model import ARIMA
from sklearn.metrics import mean_squared_error
def parser(x):
	return datetime.strptime('190'+x, '%Y-%m')


microclim_hourly=pd.read_csv('Microclimate_Hourly_threshold_0_0.csv')
par_logical=microclim_hourly.columns.str.startswith('PAR')
par=microclim_hourly.loc[:,par_logical]
row1=par.iloc[0,:]



from pandas.plotting import autocorrelation_plot
autocorrelation_plot(row1)
pyplot.show()
#1500? The first 1500 observations seems to be strongly autocorrelated... first two months

from statsmodels.tsa.arima_model import ARIMA

# fit model
model = ARIMA(row1, order=(24,2,0)) #p = the number of lag observations included in the model, d = degree of differencing, q= size of the moving average window
model_fit = model.fit(disp=0)
print(model_fit.summary())
# plot residual errors
residuals = pd.DataFrame(model_fit.resid)
residuals.plot()
pyplot.show()
residuals.plot(kind='kde')
pyplot.show()
print(residuals.describe())

'''search arima hyperparameters with python'''

# evaluate an ARIMA model for a given order (p,d,q)
def evaluate_arima_model(X, arima_order):
	# prepare training dataset
	train_size = int(len(X) * 0.66)
	train, test = X[0:train_size], X[train_size:]
	history = [x for x in train]
	# make predictions
	predictions = list()
	for t in range(len(test)):
		model = ARIMA(history, order=arima_order)
		model_fit = model.fit(disp=0)
		yhat = model_fit.forecast()[0]
		predictions.append(yhat)
		history.append(test[t])
	# calculate out of sample error
	error = mean_squared_error(test, predictions)
	return error

import warnings
warnings.filterwarnings("ignore")

# evaluate combinations of p, d and q values for an ARIMA model
def evaluate_models(dataset, p_values, d_values, q_values):
	dataset = dataset.astype('float32')
	best_score, best_cfg = float("inf"), None
	for p in p_values:
		for d in d_values:
			for q in q_values:
				order = (p,d,q)
				try:
					mse = evaluate_arima_model(dataset, order)
					if mse < best_score:
						best_score, best_cfg = mse, order
					print('ARIMA%s MSE=%.3f' % (order,mse))
				except:
					continue
	print('Best ARIMA%s MSE=%.3f' % (best_cfg, best_score))

#load dataset
microclim_hourly=pd.read_csv('Microclimate_Hourly_threshold_0_0.csv')
par_logical=microclim_hourly.columns.str.startswith('PAR')
par=microclim_hourly.loc[:,par_logical]
row1=par.iloc[0,:]

#specifiy list of parameters to test
p_values = [0, 1, 2, 4, 6, 8, 10,24,168,720]
d_values = range(0, 10)
q_values = range(0, 10)

evaluate_models(row1.values, p_values, d_values, q_values)