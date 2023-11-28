import numpy as np
import scipy
import scipy.optimize
import scipy.interpolate
import scipy.integrate
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import *
from numpy import sqrt, sin, cos, pi
from scipy.optimize import curve_fit
import sys
import os
import radialProfile

##################Function definition############################
def compute_Ra(height):
    Ra = 0
    for i in range(height.shape[0]):
        for j in range(height.shape[1]):
            Ra += abs(height[i,j])
    Ra = Ra/(nx*ny)
    return Ra

def compute_Rrms(height):
    sigma_rms_total = 0
    for i in range(height.shape[0]):
        for j in range(height.shape[1]):
            sigma_rms_total += height[i,j]**2
    sigma_rms_total = sqrt(sigma_rms_total/(nx*ny))
    return sigma_rms_total

def make_a_zero_mean_field(height):
    mean = 0
    for i in range(height.shape[0]):
        for j in range(height.shape[1]):
            mean += height[i,j]
    mean = mean/(nx*ny)
    for i in range(height.shape[0]):
        for j in range(height.shape[1]):
            height[i,j] -= mean
    return height

def plot_height_field(PATH, height, dl, SAMPLE, show=True, save=True):
    xmin = 0
    xmax = len(x)*dl
    ymin = 0
    ymax = len(y)*dl
    plt.close()
    plt.clf()
    imshow(height, interpolation=None, extent=[xmin,xmax,ymin,ymax])
    plt.xlabel("length x [" + r'$\mu$' + "m]",weight="bold")
    plt.ylabel("length y [" + r'$\mu$' + "m]", weight="bold")
    plt.title('Height field of ' + SAMPLE, weight="bold")
    cbar = plt.colorbar()
    cbar.set_label('Height [' + r'$\mu$' + 'm]', weight="bold")
    if save == True:
        plt.savefig(PATH + SAMPLE + '_Height_Field.png', dpi=400)
    if show == True:
        plt.show()

def compute_slope(height,dl):
    ny, nx = height.shape
    slope = np.zeros((ny-2,nx-2))
    for i in range(1,ny-1):
        for j in range(1,nx-1):
            d_dx = (height[i+1,j] - height[i-1,j])/(2.0*dl)
            d_dy = (height[i,j+1] - height[i,j-1])/(2.0*dl)
            slope[i-1,j-1] = sqrt(d_dx*d_dx + d_dy*d_dy)
    return slope

def plot_slope_field(PATH, slope, dl, SAMPLE, show=True,save=True):
    xmin = 0
    xmax = len(x)*dl
    ymin = 0
    ymax = len(y)*dl
    plt.close()
    plt.clf()
    imshow(slope, interpolation=None, extent=[xmin+dl,xmax-dl,ymin+dl,ymax-dl])
    plt.xlabel("length x [" + r'$\mu$' + "m]", weight="bold")
    plt.ylabel("length y [" + r'$\mu$' + "m]", weight="bold")
    plt.title('Slopes of ' + SAMPLE, weight="bold")
    cbar = plt.colorbar()
    cbar.set_label('Slope [-]', weight="bold")
    if save == True:
        plt.savefig(PATH + SAMPLE + '_Slopes_Field.png', dpi=400)
    if show == True:
        plt.show()

def read_data(FILE, PLACE):
    if PLACE == 'BASEL':
        height = np.genfromtxt(filename, skip_header=15, skip_footer=0, delimiter='";"', dtype = 'str')
    elif PLACE == 'BASEL_v':
        height = np.genfromtxt(filename, skip_header=15, skip_footer=0, delimiter='","', dtype = 'str')
    elif PLACE == 'BASEL_afm':
        height = np.loadtxt(filename, skiprows=4)
    elif PLACE == 'CEA':
        height = np.genfromtxt(filename, skip_header=15, skip_footer=0, delimiter=';', dtype = 'str')
    
    if PLACE == 'BASEL':
        file=open(filename)
        lines=file.readlines()
        value=lines[6].split(';"')[1]
        value = value.replace('"','')
        value = value.replace(',','.')
        xy_step=np.double(value)
    elif PLACE == 'BASEL_v':
        file=open(filename)
        lines=file.readlines()
        value=lines[6].split(',"')[1]
        value = value.replace('"','')
        value = value.replace(',','.')
        xy_step=np.double(value)
    elif PLACE == 'BASEL_afm':
        ny, nx = height.shape
        file=open(filename)
        lines=file.readlines()
        value=lines[1].split(': ')[1]
        xy_step=np.double(value)/nx
    elif PLACE == 'CEA':
        file=open(filename)
        lines=file.readlines()
        value=lines[6].split(';')[1]
        value = value.replace('"','')
        value = value.replace(',','.')
        xy_step=np.double(value)

    if PLACE != 'BASEL_afm':
        for i in range(height.shape[0]):
            for j in range(height.shape[1]):
                height[i,j] = height[i,j].replace('"','')
                height[i,j] = height[i,j].replace(',','.')
        height = height.astype(float)

    if PLACE == 'BASEL_afm':
        for i in range(height.shape[0]):
            for j in range(height.shape[1]):
                height[i,j] = height[i,j]*1e6

    return height, xy_step

def count_occurrences(arr, ranges):
    counts = {}
    for start, end in ranges:
        count = 0
        for sublist in arr:
            for num in sublist:
                if start <= num <= end:
                    count += 1
        counts[(start, end)] = count
    return counts   

def fractal(f, K, s):
    return K/(f**s)

def psd_func(f, A, B, C):
    return A/((1+B*B*f*f)**(C/2))

##################End of function definition############################




#def psd_func(f, psd0, xi, H):
#    return psd0/(1+abs(2*pi*f*xi)**(2*H+1))



save_picture = True
show = False


#HERE IS THE PLACE TO INPUT THE FOLDER WITH THE RESULTS AND THE PATH TO WHERE THE DATA IS STORED
results_prefix = './Results_Analysis2D/' #The "." before the / means that the path starts where this folder is
path_to_data_prefix = './Data/' #Path to where the data is stored
MAGNIFICATION = 'x50' #magnification used for the calculations
sample_list = ['039-SS', '153-SS', '673-SS', '687-SS', '689-SS', '694-SS', '737-SS','755-SS', 'Au-Si-01', 'Au-Si-02']
PLACE = 'BASEL_v'
#### SAMPLES M100, M101, M102, M103 TOMÁS #####

for sample in sample_list:

    print('Analyzing sample:', sample)
    PLACE = 'BASEL_v'
    path = path_to_data_prefix + sample + '/'
    results = results_prefix + sample + '/'
    SAMPLE = sample
    filename = path + sample + '_topography_confocal.csv'
    label = 'confocal_' + sample + MAGNIFICATION



    #Check if the results path exists
    if not os.path.exists(results):
        os.makedirs(results)


    #####################Start of the calculations#######################################################


    height, xy_step = read_data(filename, PLACE)

    print(height.shape)

    ny, nx = height.shape

    print("nx = ", nx)
    print("ny = ", ny)

    nx = height.shape[1]
    ny = height.shape[0]

    print("nx = ", nx)
    print("ny = ", ny)

    x =[0]*nx
    y =[0]*ny

    Lx = len(x)*xy_step
    Ly = len(y)*xy_step
    Lmax = max(Lx,Ly)
    dl = xy_step

    print("xy_step = ", xy_step)

    print("Lx = ", Lx)
    print("Ly = ", Ly)
    print("Lmax = ", Lmax)

    for i in range(len(x)):
        x[i]=i*xy_step

    for j in range(len(y)):
        y[j]=j*xy_step

    xmin = min(x)
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)

    height = make_a_zero_mean_field(height)

    print('min : ' + str(height.min()))
    print('max : ' + str(height.max()))

    Ra = compute_Ra(height)
    sigma_rms_total = compute_Rrms(height)

    print('sigma_rms_total = ', sigma_rms_total)
    print('Ra = ', Ra)

    plot_height_field(results,height,dl,SAMPLE,show=show,save=save_picture)
    slope = compute_slope(height,dl)
    plot_slope_field(results,slope,dl,SAMPLE,show=show,save=save_picture)



    deltaa = 0.0
    deltaq = 0.0

    for i in range(ny-2):
        for j in range(nx-2):
            deltaa += slope[i,j]

    deltaa /= (nx-2)*(ny-2)


    for i in range(ny-2):
        for j in range(nx-2):

            #deltaq += (slope[i,j] - deltaa)*(slope[i,j] - deltaa)
            deltaq += slope[i,j]*slope[i,j]

    deltaq = sqrt(deltaq/((nx-2)*(ny-2)))


    print('deltaa = ', deltaa)
    print('deltaq = ', deltaq)






    delta = arctan(slope)

    plt.close()
    plt.clf()

    imshow(delta, interpolation=None, extent=[xmin+dl,xmax-dl,ymin+dl,ymax-dl])

    plt.xlabel("length x [" + r'$\mu$' + "m]",weight="bold")
    plt.ylabel("length y [" + r'$\mu$' + "m]", weight="bold")
    plt.title('delta of ' + SAMPLE, weight="bold")

    cbar = plt.colorbar()
    cbar.set_label('Angle [rad]', weight="bold")

    if save_picture == True:
        plt.savefig(results + SAMPLE + '_Angle_Field_rad.png', dpi=400)
    if show == True:
        plt.show()


    delta_mean = np.mean(delta)





    # Exemple d'utilisation
    value_ranges = [(i, i+1) for i in range(0, 89)]

    frequency = [0]*90
    idx = 0

    occurrences = count_occurrences(delta*180/pi, value_ranges)
    for (start, end), count in occurrences.items():
        #print(f"Occurrences entre {start} et {end} : {count}")
        frequency[idx] = 100*count/((nx-2)*(ny-2))
        idx = idx + 1


    #A = np.concatenate((xpsd1D,psd1D))
    #A = transpose(A.reshape((2,len(xpsd1D))))
    np.savetxt(results + SAMPLE + '_Angle_Deg_Frequency.txt', frequency, fmt="%s")




    plt.close()
    plt.clf()

    plt.plot(frequency)

    plt.xlabel('Angle [°]', weight="bold")
    plt.ylabel('Rate of occurence [%]', weight="bold")
    plt.title('Frequency of each angle ' + SAMPLE, weight="bold")

    if save_picture == True:
        plt.savefig(results + SAMPLE + '_Angle_Deg_Frequency.png', dpi=400)
    if show == True:
        plt.show()



    plt.close()
    plt.clf()

    imshow(delta*180/pi, interpolation=None, extent=[xmin+dl,xmax-dl,ymin+dl,ymax-dl])

    plt.xlabel("length x [" + r'$\mu$' + "m]",weight="bold")
    plt.ylabel("length y [" + r'$\mu$' + "m]", weight="bold")
    plt.title('delta of ' + SAMPLE, weight="bold")

    cbar = plt.colorbar()
    cbar.set_label('Angle [°]', weight="bold")

    if save_picture == True:
        plt.savefig(results + SAMPLE + '_Angle_Field_deg.png', dpi=400)
    if show == True:
        plt.show()



    plt.close()
    plt.clf()

    imshow(delta/(pi/2), interpolation=None, extent=[xmin+dl,xmax-dl,ymin+dl,ymax-dl])

    plt.xlabel("length x [" + r'$\mu$' + "m]",weight="bold")
    plt.ylabel("length y [" + r'$\mu$' + "m]", weight="bold")
    plt.title('delta of ' + SAMPLE, weight="bold")

    cbar = plt.colorbar()
    cbar.set_label('Angle normalized [-]', weight="bold")

    if save_picture == True:
        plt.savefig(results + SAMPLE + '_Angle_Field_normalized.png', dpi=400)
    if show == True:
        plt.show()





    ft = np.fft.ifftshift(height)
    ft = np.fft.fft2(ft)
    ft = np.fft.fftshift(ft)

    #psd = (1/(Lx*Ly))*abs(ft)**2
    psd = (1/(nx*ny))*abs(ft)**2

    #psd_min = psd.min()
    #psd_max = psd.max()


    psd_min = 1e9
    psd_max = -1e9

    for i in range(psd.shape[0]):
        for j in range(psd.shape[1]):
            
            if psd[i,j] > psd_max:
                psd_max = psd[i,j]
                psd_max_i = i
                psd_max_j = j

            if psd[i,j] < psd_min:
                psd_min = psd[i,j]
                psd_min_i = i
                psd_min_j = j

    psd[psd_min_i,psd_min_j] = psd_max
            

    #print("psd_min = ", psd_min)
    #print("psd_min_i = ", psd_min_i)
    #print("psd_min_j = ", psd_min_j)
    #print("psd_max = ", psd_max)
    #print("psd_max_i = ", psd_max_i)
    #print("psd_max_j = ", psd_max_j)



    fy = [0]*height.shape[0]
    fx = [0]*height.shape[1]

    #print(len(x))
    #print(len(y))


    fx_min=0
    fx_max=1/xy_step
    fx_stp=2*(fx_max - fx_min)/(len(fx)-1)

    fy_min=0
    fy_max=1/xy_step
    fy_stp=2*(fy_max - fy_min)/(len(fy)-1)


    for i in range(len(fx)):
        #fx[i]=1/2*len(fx)/xy_step-i/xy_step
        fx[i]=-(fx_max - fx_min) + fx_min+i*fx_stp

    for j in range(len(fy)):
        #fy[j]=1/2*len(fy)/xy_step-j/xy_step
        fy[j]=-(fy_max - fy_min) + fy_min+j*fy_stp


    fxmin = min(fx)
    fxmax = max(fx)
    fymin = min(fy)
    fymax = max(fy)

    fmax = sqrt(fxmax**2 + fymax**2)



    #print('fxmin = ', fxmin)
    #print('fxmax = ', fxmax)
    #print('fymin = ', fymin)
    #print('fymax = ', fymax)
    #print('fmax = ', fmax)


    sigma_rms_total = 0

    for i in range(psd.shape[0]):
        for j in range(psd.shape[1]):
            sigma_rms_total += psd[i,j]

    sigma_rms_total = sqrt(sigma_rms_total/(nx*ny))

    print('sigma_rms_total (PSD) = ', sigma_rms_total)


    # Fit between 400 and 800nm
    f_relevant_min = 1/0.800
    f_relevant_max = 1/0.400

    sigma_rms_rel_400_800 = 0
    N = 0

    rmax = 0


    for i in range(psd.shape[0]):
        for j in range(psd.shape[1]):
            r = sqrt(fx[j]**2 + fy[i]**2)

            if r > rmax:
                rmax = r

            if r > f_relevant_min and r < f_relevant_max:
                sigma_rms_rel_400_800 += psd[i,j]
                N += 1

    sigma_rms_rel_400_800 = sqrt(sigma_rms_rel_400_800/N)

    print('sigma_rms_rel [400,800] = ', sigma_rms_rel_400_800)



    # Fit between 250 and 2500nm
    f_relevant_min = 1/2.500
    f_relevant_max = 1/0.250

    sigma_rms_rel_250_2500 = 0
    N = 0

    rmax = 0


    for i in range(psd.shape[0]):
        for j in range(psd.shape[1]):
            r = sqrt(fx[j]**2 + fy[i]**2)

            if r > rmax:
                rmax = r

            if r > f_relevant_min and r < f_relevant_max:
                sigma_rms_rel_250_2500 += psd[i,j]
                N += 1

    sigma_rms_rel_250_2500 = sqrt(sigma_rms_rel_250_2500/N)

    print('sigma_rms_rel [250,2500] = ', sigma_rms_rel_250_2500)






    #print('rmax = ', rmax)





    '''
    subdiv = 10

    radius = np.linspace(fmax/subdiv, fmax, subdiv)
    wavelength = [0]*len(radius)
    roughness_rel = [0]*len(radius)



    for k in range(len(radius)):
        wavelength[k] = 1/radius[k]
        N = 0
        for i in range(psd.shape[0]):
            for j in range(psd.shape[1]):
                r = sqrt(fx[j]**2 + fy[i]**2)

                if r < radius[k]:
                    roughness_rel[k] += psd[i,j]
                    N += 1

        roughness_rel[k] = sqrt(roughness_rel[k]/N)

    wavelength.reverse()

    #print('radius = ', radius)
    #print('wavelength = ', wavelength)

    lambda_min = min(wavelength)
    lambda_max = max(wavelength)

    print('lambda_min = ', lambda_min)
    print('lambda_max = ', lambda_max)

    plt.close()
    plt.clf()
    plt.xlabel(r'$\lambda$' + " (" + r'$\mu$' + "m" + ")",weight="bold")
    plt.ylabel(r'$\sigma_{rel}$' + " (" + r'$\mu$' + "m" + ")", weight="bold")
    plt.title('Relevant rms roughness of ' + SAMPLE, weight="bold")
    plt.plot(wavelength, roughness_rel)

    if save_picture == True:
        plt.savefig(results + SAMPLE + '_Rougness_rms_relevant.png', dpi=400)
    if show == True:
        plt.show()

    '''


    plt.close()
    plt.clf()


    #imshow(psd, interpolation=None, extent=[fxmin,fxmax,fymin,fymax])
    #imshow(psd, interpolation=None, norm=LogNorm(vmin=psd_min, vmax=psd_max), extent=[fxmin,fxmax,fymin,fymax])
    imshow(psd, interpolation=None, norm=LogNorm(), extent=[fxmin,fxmax,fymin,fymax])

    # M102
    #imshow(psd, interpolation=None, norm=LogNorm(vmin=1e-6, vmax=psd_max), extent=[fxmin,fxmax,fymin,fymax])

    plt.xlabel("frequency fx [" + r'$\mu$' + "m" + r'${}^{-1}$' + "]",weight="bold")
    plt.ylabel("frequency fy [" + r'$\mu$' + "m" + r'${}^{-1}$' + "]", weight="bold")
    plt.title('PSD of ' + SAMPLE, weight="bold")

    cbar = plt.colorbar()
    cbar.set_label('PSD [' + r'$\mu$' + 'm' + r'${}^{3}$' + ']', weight="bold")

    if save_picture == True:
        plt.savefig(results + SAMPLE + '_PSD.png', dpi=400)
    if show == True:
        plt.show()





    #psd1D = radialProfile.azimuthalAverage(psd)

    psd1D = radialProfile.radial_profile(psd)

    xpsd1D = [0]*len(psd1D)

    xpsd1D_min=1/(len(psd1D)*xy_step)
    #xpsd1D_max=1/xy_step
    xpsd1D_max=sqrt(fxmax**2 + fymax**2)
    xpsd1D_stp=(xpsd1D_max - xpsd1D_min)/(len(psd1D)-1)



    for i in range(len(xpsd1D)):
        #xpsd1D[i]=i/(xy_step*len(xpsd1D))
        #xpsd1D[i]=i/xy_step
        xpsd1D[i]=xpsd1D_min+i*xpsd1D_stp

    # unit : units^2 / frequency


    #print("xpsd1D_min = ", xpsd1D_min)
    #print("xpsd1D_max = ", xpsd1D_max)
    #print("xpsd1D_stp = ", xpsd1D_stp)




    '''
    #fmin = 1/5
    #fmax = 1/4

    #fmin = 0
    #fmax = 1/xy_step

    fmin = xpsd1D_min
    fmax = xpsd1D_max

    i_min = 0
    i_max = 0

    for i in range(len(xpsd1D)):
        if xpsd1D[i] >= fmin:
            i_min = i
            break

    for i in range(len(xpsd1D)):
        if xpsd1D[i] > fmax:
            i_max = i-1
            break

    if i_max == 0:
        i_max = len(xpsd1D)
        
    print("i_min = ", i_min)
    print("i_max = ", i_max)



    psd1D_function = scipy.interpolate.interp1d(xpsd1D,psd1D)
    res, err = scipy.integrate.quad(psd1D_function, fmin+1e-6, fmax-1e-6, limit=len(xpsd1D))

    print("res = ", res)
    print("err = ", err)

    sigma_rms_rel = 0

    i_min = 0
    i_max = len(xpsd1D)-1

    for i in range(i_min, i_max):
        sigma_rms_rel += psd1D[i]*xpsd1D[i]

    sigma_rms_rel *= (xpsd1D[i_max]-xpsd1D[i_min])/(i_max-i_min)
    sigma_rms_rel *= 2*pi
    sigma_rms_rel = sqrt(sigma_rms_rel)

    print("sigma_rms_rel = ", sigma_rms_rel)
    '''


    #popt, pcov = curve_fit(psd_func, xpsd1D, psd1D)
    '''
    popt, pcov = curve_fit(psd_func, xpsd1D, psd1D, p0=[10, 0.01, 1], bounds=(0, inf))

    psd0 = popt[0]
    xi = popt[1]
    H = popt[2]
    '''
    psd0 = 120
    xi = 1/(2*pi*30)
    H = 1

    x_fit = np.linspace(min(xpsd1D), max(xpsd1D), 1000)
    y_fit = psd_func(x_fit, psd0, xi, H)

    #print("psd0 = ", psd0)
    #print("xi = ", xi)
    #print("H = ", H)








    # Fit between 400 and 800nm
    f_relevant_min = 1/0.800
    f_relevant_max = 1/0.400

    i = 0
    i_min_xpsd1D = 0
    i_max_xpsd1D = 0

    while xpsd1D[i] < f_relevant_min:
            i_min_xpsd1D = i
            i += 1

    while i != len(xpsd1D) and xpsd1D[i] < f_relevant_max:
            i_max_xpsd1D = i
            i += 1


    #print("i_min_xpsd1D = ", i_min_xpsd1D)
    #print("i_max_xpsd1D = ", i_max_xpsd1D)


    popt, pcov = curve_fit(fractal, xpsd1D[i_min_xpsd1D:i_max_xpsd1D], psd1D[i_min_xpsd1D:i_max_xpsd1D])

    K_400_800 = popt[0]
    s_400_800 = popt[1]

    print("K_400_800 = ", K_400_800)
    print("s_400_800 = ", s_400_800)





    # Fit between 250 and 2500nm
    f_relevant_min = 1/2.500
    f_relevant_max = 1/0.250

    i = 0
    i_min_xpsd1D = 0
    i_max_xpsd1D = 0

    while xpsd1D[i] < f_relevant_min:
            i_min_xpsd1D = i
            i += 1

    while i != len(xpsd1D) and xpsd1D[i] < f_relevant_max:
            i_max_xpsd1D = i
            i += 1

    popt, pcov = curve_fit(fractal, xpsd1D[i_min_xpsd1D:i_max_xpsd1D], psd1D[i_min_xpsd1D:i_max_xpsd1D])

    K_250_2500 = popt[0]
    s_250_2500 = popt[1]

    print("K_250_2500 = ", K_250_2500)
    print("s_250_2500 = ", s_250_2500)








    plt.close()
    plt.clf()
    plt.xlabel("frequency [" + r'$\mu$' + "m" + r'${}^{-1}$' + "]",weight="bold")
    plt.ylabel("PSD [" + r'$\mu$' + "m" + r'${}^3$' + "]", weight="bold")
    plt.title('Angularly avergaged PSD of ' + SAMPLE, weight="bold")
    plt.xscale("log")
    plt.yscale("log")
    plt.plot(xpsd1D, psd1D, label=label)#, linestyle='-', marker='o')
    #plt.plot(x_fit, y_fit)

    plt.legend()

    from matplotlib.patches import Rectangle

    #left, bottom, width, height = (xpsd1D[i_min_xpsd1D], min(psd1D), xpsd1D[i_max_xpsd1D]-xpsd1D[i_min_xpsd1D], max(psd1D)-min(psd1D))

    bottom, top = plt.gca().get_ylim()
    left, bottom, width, height = (f_relevant_min, bottom, f_relevant_max-f_relevant_min, top-bottom)

    rect=Rectangle((left,bottom),width,height, alpha=0.1, facecolor="red")
    plt.gca().add_patch(rect)


    if save_picture == True:
        plt.savefig(results + SAMPLE + '_PSD1D.png', dpi=400)
    if show == True:
        plt.show()


    A = np.concatenate((xpsd1D,psd1D))
    A = transpose(A.reshape((2,len(xpsd1D))))
    np.savetxt(results + SAMPLE + '_psd1D_data.txt', A, fmt="%s")



    with open(results + SAMPLE + '_res.txt', 'w') as file:
        file.write('sigma_rel 400-800nm\t' + str(sigma_rms_rel_400_800) + '\n')
        file.write('sigma_rel 250-2500nm\t' + str(sigma_rms_rel_250_2500) + '\n')
        file.write('Ra\t' + str(Ra) + '\n')
        file.write('Rs\t' + str(deltaa) + '\n')
        file.write('RMS height\t' + str(sigma_rms_total) + '\n')
        file.write('RMS slopes\t' + str(deltaq) + '\n')
        file.write('K 400-800nm\t' + str(K_400_800) + '\n')
        file.write('s 400-800nm\t' + str(s_400_800) + '\n')
        file.write('K 250-2500nm\t' + str(K_250_2500) + '\n')
        file.write('s 250-2500nm\t' + str(s_250_2500) + '\n')
        file.write('delta_m\t' + str(delta_mean) + '\n')






