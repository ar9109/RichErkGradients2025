function [myoutput] = odr_scipy(x,y)
%ODR_SCIPY Summary of this function goes here
%   Detailed explanation goes here
odr = py.importlib.import_module('scipy.odr');
np = py.importlib.import_module('numpy');
if count(py.sys.path,'D:\BIO\PhD\ditalia\zebrafish\github\ditalia-zebrafish\code\functions\py') == 0
    insert(py.sys.path,int32(0),'D:\BIO\PhD\ditalia\zebrafish\github\ditalia-zebrafish\code\functions\py');
end
mymod = py.importlib.import_module('mymod');
py.importlib.reload(mymod);

% construct fit
xdata = np.asarray(x,'float64');
ydata = np.asarray(y,'float64');

f = mymod.fit_func_py;
power_py = odr.Model(mymod.fit_func_py);
mydata = odr.Data(xdata, ydata);
myodr = odr.ODR(mydata, power_py, np.array([1, 1]));
myoutput = myodr.run();

end

