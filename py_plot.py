import matplotlib.pyplot as plt
import numpy as np
import statistics as st
import scipy.optimize as opt
import scipy
import math

accuracy = 1000
DESCEND_ACCEL = 0.2
SYMBOLS_AFTER_POINT = 3
DIRECTION_STEPS = 10
DESCENT_RATE_MIN = 1e-20
DEBUG = False

I0 = 200

def get_data(num):
    return float(num)
def get_y_data(num):
    print(num)
    return float(num)#-math.log(float(num)/I0)
def get_yerr(num):
    return float(num)

def find_k_and_b(x, y):
    k = (st.mean(x*y)-st.mean(x)*st.mean(y))/(st.mean(x*x) - st.mean(x)**2)
    return k, (st.mean(y)-k*st.mean(x))
def find_k(x, y):
    return st.mean(x*y)/st.mean(x*x)

#params[0] - k
#params[1] - b
def parabol_approx(x, a, b, c):
    return a*x**2+b*x+c

def linear_approx(x, k, b):
    return k*x+b

def const_approx(x, b):
    return b

def gauss(x, A, x0, sigma): 
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    
# def const_approx(x, params):
#     if len(params) == 1:
#         return params[0]

def chi_square(y_pred, y_data, yerr):
    # print("---------Chi_square debug---------")
    # print("Prediction: {}".format(y_pred))
    # print("Data: {}".format(y_data))
    # print("Error: {}".format(yerr))
    # print("Delta: {}".format(((np.array(y_pred)-np.array(y_data))/(np.array(yerr)))**2))
    # print("===============END================")
    return np.sum( ((np.array(y_pred)-np.array(y_data))/(np.array(yerr)))**2 )

# def approx_chi_sqaure(x_data, y_data, yerr, params=[1, 0], y_func=linear_approx, res_func=chi_square, descent_rate=[0.000001, 0.00001], iteration_num=120000, delta_param=1e-9):
#     y_pred = []
#     gradient = [param for param in params]
#     direction = [0 for _ in params]
#     direction_time = [0 for _ in params]
#     res = -1
#     for i in range(iteration_num):
#         y_pred = [y_func(x_data[k], params) for k in range(len(y_data))]
#         res = res_func(y_pred, y_data, yerr)
#         #plt.errorbar(x_data, y_data, yerr = yerr, fmt = '.')
#         #plt.errorbar(x_data, y_pred, fmt = '.')
#         #plt.show()

#         if(DEBUG):
#             print("\n----------------\nIteration {}, res = {}".format(i, res))

#         for j in range(len(params)):
#             temp_params = [param for param in params]
#             if(DEBUG):
#                 print("Params: {} vs Copy: {}".format(params, temp_params))
#             temp_params[j] += delta_param
#             temp_res = res_func( [y_func(x_data[k], temp_params) for k in range(len(y_data))], y_data, yerr )
#             if(DEBUG):
#                 print("Res is {} for params: {}".format(temp_res, temp_params))
#                 print("Delta res is {}, gradient[{}] = {}".format(temp_res-res, j, (temp_res - res)/delta_param))
#             gradient[j] = (temp_res - res) / delta_param
#             if(direction[j]*(temp_res - res) < 0):
#                 descent_rate[j] = max(descent_rate[j]*DESCEND_ACCEL, DESCENT_RATE_MIN)
#                 direction_time[j] = 0
#             elif direction_time[j] >= DIRECTION_STEPS:
#                 direction_time[j] = 0
#                 descent_rate[j] = min(descent_rate[j]/DESCEND_ACCEL, 1)
#             else:
#                 direction_time[j]+=1
#             direction[j] = 1 if (temp_res - res) >= 0 else -1
        
#         for j in range(len(params)):
#             params[j] -= descent_rate[j]*gradient[j]
        
#     params_err = [param/accuracy/100 for param in params]
#     for j in range(len(params)):
#         temp_params = [param for param in params]
#         for k in range(0, accuracy*100):
#             temp_res = res_func( [y_func(x_data[g], temp_params) for g in range(len(y_data))], y_data, yerr )
#             print("Counting sigma for param[{}]: k = {}/{}\ntemp_res = {}\nres = {}".format(j, k, accuracy*100, temp_res, res))
#             if((temp_res - res)**2 > 1):
#                 params_err[j] *= k
#                 break
#             temp_params[j] += params_err[j]

#     return params, res, params_err

names = []
def make_format_array(popt, perr, i=math.nan):
    m_popt = np.array(popt)
    m_perr = np.array(perr)
    format_arr = []
    if not math.isnan(i):
        format_arr += [names[int(i)]]
    for k in range(len(m_popt)):
        mnoj = 0
        while(math.fabs(m_popt[k]) < 0.1):
            m_popt[k]*=10
            m_perr[k]*=10
            mnoj-=1
        while(math.fabs(m_popt[k]) > 10):
            m_popt[k]/=10
            m_perr[k]/=10
            mnoj+=1
        if(mnoj == 0):
            format_arr += [round(m_popt[k], SYMBOLS_AFTER_POINT)]
            format_arr += [round(m_perr[k], SYMBOLS_AFTER_POINT)]
        else:
            format_arr += ["{} e{}".format(round(m_popt[k], SYMBOLS_AFTER_POINT), mnoj)]
            format_arr += ["{} e{}".format(round(m_perr[k], SYMBOLS_AFTER_POINT), mnoj)]
    return format_arr

file = open("data.csv")
data = []
for line in file:
    data.append(''.join(line.replace('\t', ',').replace('\n', ',\n').split()).split(","))

# for i in range(len(data)):
#     for j in range(len(data[i])):
#         if(data[i][j] != ''):
#             data[i][j] = float(data[i][j])
# data.sort()
#for dat in data:
 #    print(dat)

k = 1
b = 0
y_data = 1
x_data = 1
yerr = 1
xerr = 1
MNK = False
CHI_SQUARE = True
CHI_CONST =  False
CHI_SQUARE_PARAB = False
CHI_SQUARE_GAUSS = False
AKIMA_INTERPOL = False
DRAW_ALL = False
MNK_ZERO_IS_ZERO = False
x_range = (0, 100)

SUM = 0


for i in range(0, 2, 1):

    y_data = []
    yerr = []

    data_len = len(data)
    for j in range(0, data_len):
        if data[j][i] == '':
            data_len = j
            break

    for y in range(0, data_len):
            y_data += [get_y_data(data[y][i+1])]
            yerr += []
            # y_data += [math.sqrt(get_y_data(data[y][i+1]))/math.sqrt(get_data(data[y][i]))**3*10**6]
            # yerr += [get_yerr(data[y][i+2])/2/math.sqrt(get_y_data(data[y][i+1]))/math.sqrt(get_data(data[y][i]))**3*10**6]
    x_data = []
    xerr = []
    for x in range(0, data_len):
            x_data += [get_data(data[x][0])]
            xerr += [0.1]
            # x_data += [math.sqrt(get_data(data[x][i])**2 + mc_square**2) - mc_square]
            # xerr += [x_data[x-1]*((0.01*K/get_data(data[x][i]))+(3/K))]

    # if(MNK_ZERO_IS_ZERO):
    #      k = find_k(np.array(x_data), np.array(y_data))
    # else:
    #     k, b = find_k_and_b(np.array(x_data), np.array(y_data))
    #yerr = [math.log(1.03) for y in range(0, len(y_data), 1)]
    # xerr = [0.01*K for x in range(0, len(x_data))]

    if AKIMA_INTERPOL:
        offset = 1
        offset2 = 100
        interpol = scipy.interpolate.Akima1DInterpolator(x_data, y_data, axis=0)
        X = [x/accuracy/1 for x in range(int(min(x_data[offset:offset2])*accuracy) , int(max(x_data[offset:offset2])*accuracy), 1)]
        plt.plot(X, interpol(X))
    #print(x_data)
    plt.errorbar(x_data, y_data, yerr = yerr, xerr = xerr, fmt = '.')

    offset_arr =  [0, 0, 0, 0, 0, 0]
    offset2_arr = [8, 7, 6, 5, 5, 3]
    #lambda_arr = [5852]
    if(CHI_SQUARE):
        offset = offset_arr[int(i/2)]#7
        offset2 = offset2_arr[int(i/2)]#20
        #params, result, params_err = approx_chi_sqaure(x_data[offset2:], y_data[offset2:], yerr = yerr[offset2:], params=[4, -23.4], y_func=linear_approx, descent_rate=[1e-5, 1e-5])
        popt, pcov = opt.curve_fit(linear_approx, x_data[offset:offset2], y_data[offset:offset2], sigma=yerr[offset:offset2], absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        print(popt, perr)
        X = [x/accuracy/1 for x in range(int(min(x_data[offset:offset2])*accuracy) , int(max(x_data[offset:offset2])*accuracy), 1)]
        y = [linear_approx(x, popt[0], popt[1]) for x in X]

        #plt.errorbar(x_data[offset:offset2], y_data[offset:offset2], yerr = yerr[offset:offset2], xerr = xerr[offset:offset2], fmt = '.')
        plt.plot(X, y, label = 'k = {} ± {},  b = {} ± {}, χ^2 = {}'.format(*make_format_array(popt, perr), round(chi_square([linear_approx(x, popt[0], popt[1]) for x in x_data[offset:offset2]], y_data[offset:offset2], yerr[offset:offset2]), SYMBOLS_AFTER_POINT)))

    if(CHI_SQUARE_PARAB):
        offset = 0
        offset2 = 25
        #params, result, params_err = approx_chi_sqaure(x_data[offset2:], y_data[offset2:], yerr = yerr[offset2:], params=[4, -23.4], y_func=linear_approx, descent_rate=[1e-5, 1e-5])
        popt, pcov = opt.curve_fit(parabol_approx, x_data[offset:offset2], y_data[offset:offset2], sigma=yerr[offset:offset2], absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        print(popt, perr)
        X = [x/accuracy/1 for x in range(int(min(x_data[offset:offset2])*accuracy) , int(max(x_data[offset:offset2])*accuracy), 1)]
        y = [parabol_approx(x, popt[0], popt[1], popt[2]) for x in X]

        plt.errorbar(x_data[offset:offset2], y_data[offset:offset2], yerr = yerr[offset:offset2], xerr = xerr[offset:offset2], fmt = '.')
        plt.plot(X, y, label = 'a = {} ± {},  b = {} ± {},  c = {} ± {}, χ^2 = {}'.format(*make_format_array(popt, perr), round(chi_square([parabol_approx(x, popt[0], popt[1], popt[2]) for x in x_data[offset:offset2]], y_data[offset:offset2], yerr[offset:offset2]), SYMBOLS_AFTER_POINT)))


    if(CHI_SQUARE_GAUSS):
        offset = 25
        offset2 = 34
        #params, result, params_err = approx_chi_sqaure(x_data[offset2:], y_data[offset2:], yerr = yerr[offset2:], params=[4, -23.4], y_func=linear_approx, descent_rate=[1e-5, 1e-5])
        popt, pcov = opt.curve_fit(gauss, x_data[offset:offset2], y_data[offset:offset2], p0=[7, 1013.5, 50], sigma=yerr[offset:offset2], absolute_sigma=True)
        #popt = [7, 1013.5, 50]
        perr = np.sqrt(np.diag(pcov))
        #perr = [0, 0, 0]
        print(popt, perr)
        X = [x/accuracy/1 for x in range(int(min(x_data[offset:offset2])*accuracy) , int(max(x_data[offset:offset2])*accuracy), 1)]
        y = [gauss(x, popt[0], popt[1], popt[2]) for x in X]

        plt.errorbar(x_data[offset:offset2], y_data[offset:offset2], yerr = yerr[offset:offset2], xerr = xerr[offset:offset2], fmt = '.')
        plt.plot(X, y, label = 'a = {} ± {},  b = {} ± {},  c = {} ± {}, χ^2 = {}'.format(*make_format_array(popt, perr), round(chi_square([gauss(x, popt[0], popt[1], popt[2]) for x in x_data[offset:offset2]], y_data[offset:offset2], yerr[offset:offset2]), SYMBOLS_AFTER_POINT)))


    if(CHI_CONST):
        offset = 0
        offset2 = 10
        popt, pcov = opt.curve_fit(const_approx, x_data[offset:offset2], y_data[offset:offset2], sigma=yerr[offset:offset2], absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        print(popt, perr)
        X = [x/accuracy/1 for x in range(int(min(x_data[offset:offset2])*accuracy) , int(max(x_data[offset:offset2])*accuracy), 1)]
        y = [const_approx(x, popt[0]) for x in X]

        plt.errorbar(x_data[offset:offset2], y_data[offset:offset2], yerr = yerr[offset:offset2], xerr = xerr[offset:offset2], fmt = '.')
        plt.plot(X, y, label = 'b = {} ± {}, χ^2 = {}'.format(*make_format_array(popt, perr), round(chi_square([const_approx(x, popt[0]) for x in x_data[offset:offset2]], y_data[offset:offset2], yerr[offset:offset2]), SYMBOLS_AFTER_POINT)))

# if(MNK):

#     size = len(y_data)
#     offset = 0
#     X = [x/accuracy/1 for x in range(x_range[0]*accuracy, x_range[1]*accuracy, 1)]
#     y = [b + k*x for x in X]

#     arr = [(-1)**(int(x/size*2)) for x in range(0, size)]
#     a = [y_data[x] + yerr[x]*arr[x - offset] for x in range(offset, offset + size)]
#     c = [y_data[x] - yerr[x]*arr[x - offset]  for x in range(offset, offset + size)]
#     d = [y_data[x] + yerr[x] for x in range(offset, offset + size)]
#     e = [y_data[x] - yerr[x] for x in range(offset, offset + size)]

#     #plt.errorbar([get_data(data[x][0]) for x in range(offset, offset + size)], a, yerr = yerr, xerr = 0, fmt = '.')

#     a1, _ = find_k_and_b(np.array([x_data[x] for x in range(offset, offset + size)]), np.array(a))
#     c1, _ = find_k_and_b(np.array([x_data[x] for x in range(offset, offset + size)]), np.array(c))
#     _, d1 = find_k_and_b(np.array([x_data[x] for x in range(offset, offset + size)]), np.array(d))
#     _, e1 = find_k_and_b(np.array([x_data[x] for x in range(offset, offset + size)]), np.array(e))
    
#     if not MNK_ZERO_IS_ZERO:
#         plt.plot(X, y, label = 'k = {} ± {},  b = {} ± {}'.format(int(k*accuracy)/accuracy, int(abs(c1 - a1)/2*accuracy)/accuracy, int(b*accuracy)/accuracy, int(abs(d1 - e1)/2*accuracy)/accuracy))
#     else:
#         plt.plot(X, y, label = 'k = {} ± {}'.format(int(k*accuracy)/accuracy, int(abs(c1 - a1)/2*accuracy)/accuracy))


plt.grid(True)
plt.ylabel("V_prob, mV")#"N^(1/2)/p^(3/2)*10^6")#
plt.xlabel("V_res, mV")
plt.title("V_prob(V_res)")
plt.legend(loc='lower right')
plt.show()
