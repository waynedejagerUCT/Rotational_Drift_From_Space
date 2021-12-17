import numpy as np
import xarray as xr
import datetime
import metpy.calc as mpcalc
from metpy.units import units
import scipy.spatial as spatial
import pandas as pd


def region_grid(product, sea):
    '''
    Create an array of the subregion over which to iterate the search algorithm
    '''
    if (product=='multi-oi') or (product=='amsr2-gw1') or (product=='ascat-metopA') or (product=='ssmi-f15') or (product=='ssmis-f17') or (product=='ssmis-f18'):
        if sea=='weddell':
            xc_max =  62.5*4
            xc_min = -62.5*40
            yc_max =  62.5*60
            yc_min =  62.5*12
        elif sea=='ross':
            xc_max = 'unknown'
            xc_min = 'unknown'
            yc_max = 'unknown'
            yc_min = 'unknown'
        elif sea=='amundsen':
            xc_max = 'unknown'
            xc_min = 'unknown'
            yc_max = 'unknown'
            yc_min = 'unknown'
        elif sea=='atlantic':
            xc_max =  62.5*34
            xc_min = -62.5*40
            yc_max =  62.5*60
            yc_min =  62.5*12
        else:
            print ('**Error finding sea**')
        
        xc = np.arange(-3875, 3912.5, 62.5)
        yc = np.arange(4250, -3912.5, -62.5)
        xc1 = np.arange(xc_min, xc_max+62.5, 62.5)
        yc1 = np.arange(yc_max, yc_min-62.5, -62.5)

    else:
        print ('**Error finding product**')
    
    return xc, yc, xc1, yc1


if __name__ == '__main__':
    '''
    EUMETSAT OSI-405-c Low Res products: multi-oi, amsr2-gw1, ascat-metopA, ssmi-f15, ssmis-f17, ssmis-f18
    '''
    product           = 'ssmis-f18'
    sea               = 'atlantic'
    radius            = 400
    isnan_threshold   = 0.1
    version           = 'v003'
    delta_1day        = datetime.timedelta(days = 1)
    
    xc, yc, xc1, yc1  = region_grid(product, sea)
    
    nan_arr           = np.empty((131,125))
    nan_arr[:]        = np.nan

    
    for year in range(2013,2021):
        print()
        print ('[' + str(year) + ']')
        print()
        
        date_0 = datetime.datetime(year, 6, 1,12,0,0)
        date_1 = datetime.datetime(year,10,31,12,0,0)
        date   = date_0
        
        df     = pd.DataFrame(columns = ['date_0','date_1','nan_percentage','x','y','mean_vort','std_vort','mean_vort_uncert','std_vort_uncert'])
        
        while date <= date_1:
            print('Date: ' + str(date))
            
            day0  = date - delta_1day - delta_1day
            day1  = date 
            
            try:
                fname     = 'ice_drift_sh_polstere-625_'+product+'_'+str(day0.year)+str(day0.month).zfill(2)+str(day0.day).zfill(2)+'1200-'+str(day1.year)+str(day1.month).zfill(2)+str(day1.day).zfill(2)+'1200.nc'
                fdir      = 'my_directory/'+product+'/'+str(day1.year)+'/'+str(day1.month).zfill(2)+'/'
                ds        = xr.open_dataset(fdir+fname)
                flag      = ds.status_flag.values[0]
                flag      = (flag>=20)
                #extract drift values: original units in kilometers per 48hrs, modified units in meters per second
                dX        = ds.dX[0].values
                dY        = ds.dY[0].values
                dX[~flag] = np.nan
                dY[~flag] = np.nan
                U         = ((dX)*1000)/(2*24*60*60)
                V         = ((dY)*1000)/(2*24*60*60)
                ICE_VORT  = mpcalc.vorticity(U*units.meter/units.second, V*units.meter/units.second, 62.5*units.kilometer, -62.5*units.kilometer)
                #extract drift uncertaity array
                if year >= 2017:
                    uncert_available                 = True
                    uncert_drift                     = ds.uncert_dX_and_dY[0].values*1000
                    uncert_drift[np.isnan(ICE_VORT)] = np.nan
                    L                                = 62.5*1000
                    delT                             = (2*24*60*60)
                    factor                           = (2)/((L**2)*(delT**2))
                    uncert_vort                      = np.sqrt(factor*np.square(uncert_drift))

                else:
                    uncert_available                 = False
                    uncert_vort                      = nan_arr
                
                points     = np.transpose([np.tile(xc, len(yc)), np.repeat(yc, len(xc))])
                point_tree = spatial.cKDTree(points)    
        
                for x in xc1:
                    for y in yc1:
            
                        points_region1  = point_tree.data[point_tree.query_ball_point([x, y], radius)]
                        VortArr   = []
                        UncertArr = []
        
                        for i in points_region1:
                            value_xc       = i[0]
                            value_yc       = i[1]
                    
                            index_xc       = np.where(xc == value_xc)[0][0]
                            index_yc       = np.where(yc == value_yc)[0][0]
                            VortArr.append(ICE_VORT[index_yc][index_xc].magnitude)
                            UncertArr.append(uncert_vort[index_yc][index_xc])
                     
                        if sum(np.isnan(VortArr))/len(VortArr) <= isnan_threshold:
                            MeanVort   = np.nanmean(VortArr)
                            std_vort   = np.nanstd(VortArr)
                            
                            if uncert_available == True:
                                MeanUncert = np.nanmean(UncertArr)
                                std_uncert = np.nanstd(UncertArr)
                            else:
                                MeanUncert = np.nan
                                std_uncert = np.nan
                            
                            df = df.append({'date_0'           : day0,
                                            'date_1'           : day1,
                                            'nan_percentage'   : (sum(np.isnan(VortArr))/len(VortArr))*100,
                                            'x'                : x,
                                            'y'                : y,
                                            'mean_vort'        : MeanVort,
                                            'std_vort'         : std_vort,
                                            'mean_vort_uncert' : MeanUncert,
                                            'std_vort_uncert'  : std_uncert}, ignore_index=True)
                            
            except(OSError):
                print('[X] OSError: File not found')
                pass
            
            print('---------------------------------')
            date = date + delta_1day
            
        #Save to csv
        df.to_csv('/my_directory/my_dataframe.csv')
