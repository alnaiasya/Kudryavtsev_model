import numpy as np
from osgeo import gdal

# функция чтения растров в формате GeoTIFF
def read_tif(path):
    file = gdal.Open(path)
    array = file.ReadAsArray()
    gtf = list(file.GetGeoTransform())
    proj = file.GetProjection()
    nodata_value = file.GetRasterBand(1).GetNoDataValue()
    return array, gtf, proj, nodata_value

# функция записи растров в формате GeoTIFF
def write_tif(path, array, gtf, proj, nodata_value):
    driver = gdal.GetDriverByName('GTiff')
    gtiff = driver.Create(path, array.shape[1], array.shape[0], 1, gdal.GDT_Float32)
    gtiff.SetGeoTransform(gtf)
    gtiff.SetProjection(proj)
    gtiff.GetRasterBand(1).WriteArray(array)
    gtiff.GetRasterBand(1).SetNoDataValue(nodata_value)
    gtiff.FlushCache()
    gtiff = None
    return

### Блок обработки данных SoilGrids

# чтение входных данных, загруженных из SoilGrids с разрешением 1000 м и осреднённых по профилю
theta_sand, gtf, proj, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\sand_v1.tif") # массовая доля песка
theta_silt, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\silt_v1.tif") # массовая доля ила
theta_clay, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\clay_v1.tif") # массовая доля глины
theta_org, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\soc_v1.tif") # Soil Organic Carbon (содержание органического углерода)
omega_water, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\vwc_v1.tif") # объёмное влагосодержание
rho_s, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\bdod_v1.tif") # плотность грунта

# замена пустых значений на NaN
theta_sand = np.where(theta_sand < 0, np.nan, theta_sand)
theta_silt = np.where(theta_silt < 0, np.nan, theta_silt)
theta_clay = np.where(theta_clay < 0, np.nan, theta_clay)
theta_org = np.where(theta_org < 0, np.nan, theta_org)
omega_water = np.where(omega_water < 0, np.nan, omega_water)
rho_s = np.where(rho_s < 0, np.nan, rho_s)

#### Табличные значения переменных (константы)

# коэффициент теплопроводности талого сухого грунта в Вт/м·°C
k_sand_t = 1.05
k_silt_t = 1.05
k_clay_t = 0.9
k_peat_t = 0.35

# коэффициент теплопроводности мёрзлого сухого грунта в Вт/м·°C
k_sand_f = 1.25
k_silt_f = 1.25
k_clay_f = 1.15
k_peat_f = 0.8

# объёмная теплоёмкость грунта в Дж/кг·°C
c_sand = 690
c_silt = 730
c_clay = 900
c_peat = 200

# плотность воды
rho_water = 1000 

#### Перевод единиц измерения и обработка выбросов

# массовая доля фракции в почве в долях единицы (сумма песка, ила и глины = 1)
theta_sand = theta_sand / 1000
theta_silt = theta_silt / 1000
theta_clay = theta_clay / 1000
theta_org = theta_org * 1.724 / 10000 # перевод SOC в SOM, значения в долях единицы, общепринятый коэфф. - 1.724

# объёмное влагосодержание в м³/м³ (исходные данные были в (значение) * 0.1 cm³/cm³)
omega_water = omega_water / 1000

# плотность почвы в кг/м³ (была в cg/cm³)
rho_s = rho_s * 10 

# обработка выбросов через условие: theta_sand + theta_silt + theta_clay == 1
theta_sum = theta_sand + theta_silt + theta_clay
theta_sum_rounded = np.round(theta_sum)
theta_sand[theta_sum_rounded != 1] = np.nan
theta_silt[theta_sum_rounded != 1] = np.nan
theta_clay[theta_sum_rounded != 1] = np.nan

#### Вычисление теплофизических характеристик грунтов

k_mineral_t = (k_sand_t ** theta_sand) * (k_silt_t ** theta_silt) * (k_clay_t ** theta_clay)
k_mineral_f = (k_sand_f ** theta_sand) * (k_silt_f ** theta_silt) * (k_clay_f ** theta_clay)

k_s_t = (k_mineral_t ** (1 - theta_org)) * (k_peat_t ** theta_org)
k_s_f = (k_mineral_f ** (1 - theta_org)) * (k_peat_f ** theta_org)

lambda_thawed_al = (k_s_t ** (1 - omega_water)) * (0.54 ** omega_water) # коэффициент теплопроводности для талых пород в Вт/м·°C
lambda_frozen_al = (k_s_f ** (1 - omega_water)) * (2.35 ** omega_water) # коэффициент теплопроводности для мёрзлых пород в Вт/м·°C


c_mineral = (c_sand * theta_sand) + (c_silt * theta_silt) + (c_clay * theta_clay)
c_s = (1 - theta_org) * c_mineral + (theta_org * c_peat)

c_vol_thawed_al = c_s * rho_s + 4190 * omega_water * (rho_water / rho_s) # объёмная теплоёмкость талых пород Дж/м³·°C
c_vol_frozen_al = c_s * rho_s + 2025 * omega_water * (rho_water / rho_s) # объёмная теплоёмкость мёрзлых пород Дж/м³·°C

# Перевод значений теплоёмкости в Вт·ч/(м³·°C)
c_vol_thawed_al = c_vol_thawed_al / 3600
c_vol_frozen_al = c_vol_frozen_al / 3600

q_phase_transition = 335200 * omega_water * rho_s # теплота фазовых переходов воды в Дж/м³
q_phase_transition = q_phase_transition / 3600 # перевод в Вт·ч/м³

lambda_pl = lambda_frozen_al # коэффициент теплопроводности для мёрзлых пород, подстилающих СТС
lambda_th = lambda_thawed_al # коэффициент теплопроводности для талых пород, подстилающих СМС
c_vol_pl = c_vol_frozen_al # объёмная теплоёмкость для мёрзлых пород, подстилающих СТС
c_vol_th = c_vol_thawed_al # объёмная теплоёмкость для талых пород, подстилающих СМС

### Чтение и обработка прочих входных данных

# показатели, связанные с t воздуха

dd_above_zero, gtf, proj, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\above_zero_sum_v1.tif") # сумма градусочасов для тёплого периода
dd_below_zero, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\below_zero_sum_v1.tif") # сумма градусочасов для холодного периода
tau_summer, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\hours_above_zero_v1.tif") # длительность тёплого периода в часах
tau_winter, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\hours_below_zero_v1.tif") # длительность холодного периода в часах
nt, _, _, _ = read_tif(r"E:\graduate_work\interpolation\Kriging_nt_c_clip.tif") # радиационная поправка для тёплого периода

omega_summer = np.where(dd_above_zero <= 0, np.nan, dd_above_zero)
omega_winter = np.where(dd_below_zero >= 0, np.nan, dd_below_zero)
tau_summer = np.where(tau_summer  <= 0, np.nan, tau_summer)
tau_winter = np.where(tau_winter <= 0, np.nan, tau_winter)

omega_summer = omega_summer * nt # сумма градусочасов на дневной поверхности для летнего периода с учётом радиационной поправки

T = 8760 * 22 # длительность периода (в часах)

# показатели, связанные с напочвенным покровом

snow_depth, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\snow_depth_v1.tif") # высота снежного покрова
snow_density, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\snow_density_v1.tif") # плотность снежного покрова
r_summer_veg, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\r_summer_glwd_v1.tif") # термическое сопротивление растительности в тёплый период
r_winter_veg, _, _, _ = read_tif(r"E:\graduate_work\v1_calc\projected_rasters\r_winter_glwd_v1.tif") # термическое сопротивление растительности в холодный период

snow_depth = np.where(snow_depth >= 1.3, np.nan, snow_depth) # ограничение высоты снега (ср. высота по метеостанциям * 2 = 1.3 м)
snow_depth = np.where(snow_depth  <= 0, np.nan, snow_depth)
snow_density = np.where(snow_density <= 0, np.nan, snow_density) 
r_summer_veg = np.where(r_summer_veg < 0, np.nan, r_summer_veg)
r_winter_veg = np.where(r_winter_veg < 0, np.nan, r_winter_veg)

lambda_snow = 0.02093 + 1.01 * (10**-3) * snow_density # вычисление коэффициента эффективной теплопроводности снега в Вт/(м³·°C) по ф-ле Проскурякова
r_snow = snow_depth / lambda_snow # среднее за зиму термическое сопротивление снежного покрова в м²·°C/Вт

### Вычисления для мёрзлых пород (СТС)

mu = 1.8 * np.sqrt(lambda_pl * c_vol_pl * tau_summer) 
s = lambda_thawed_al * r_summer_veg
alpha = lambda_thawed_al / lambda_frozen_al
beta = (omega_summer * c_vol_thawed_al) / (2 * tau_summer)
rho = r_snow + r_winter_veg - alpha * r_summer_veg # общее термическое сопротивление покровов
ground_t_zero = (alpha * omega_summer + omega_winter)/T # среднегодовая температура грунта на подошве СТС при условии отсутствия теплоизолирующих слоёв

# коэффициенты для уравнения относительно В: общей величины летнего теплооборота в породах (= зимнему значению)
a = 1 + mu * (rho / T) - 2 * s * beta * (r_summer_veg / omega_summer)
b = mu * ground_t_zero + 2 * s * (q_phase_transition + 2 * beta)
c = 2 * lambda_thawed_al * omega_summer * (q_phase_transition + beta)

B = (np.sqrt(b**2 + 4 * a * c) - b) / (2 * a) # общая величина летнего теплооборота в породах
xi = 2 * lambda_thawed_al * ((omega_summer / B) - r_summer_veg) # мощность СТС
t_xi = ground_t_zero + B * (rho / T) # среднегодовая температура пород на подошве сезонно-талого слоя

### Вычисления для талых пород (СМС)

r_winter = r_snow + r_winter_veg
mu_th = 1.8 * np.sqrt(lambda_th * c_vol_th * tau_winter)
s_th = lambda_frozen_al * r_winter
beta_th = - (omega_winter * c_vol_frozen_al) / (2 * tau_winter)
alpha_th = lambda_thawed_al / lambda_frozen_al
rho_th =  (r_winter / alpha_th) - r_summer_veg # общее термическое сопротивление покровов
ground_t_zero_th = (omega_summer + (omega_winter / alpha_th)) / T # среднегодовая температура грунта на подошве СМС при условии отсутствия теплоизолирующих слоёв

# коэффициенты для уравнения относительно В: общей величины зимнего теплооборота в породах
a_th = mu_th * (rho_th / T) - 2 * s_th * beta_th * (r_winter / omega_winter) - 1
b_th = 2 * s_th * (q_phase_transition + 2 * beta_th) - mu_th * ground_t_zero_th
c_th = -2 * lambda_frozen_al * omega_winter * (q_phase_transition + beta_th)

B_th = (b_th - np.sqrt(b_th**2 - 4 * a_th * c_th)) / (2 * a_th) # общая величина зимнего теплооборота в породах
xi_th = -2 * lambda_frozen_al * ((omega_winter / B_th) + r_winter) # мощность СМС
t_xi_th = ground_t_zero_th + B_th * rho_th / T # среднегодовая температура пород на подошве сезонно-мёрзлого слоя


### Сохранение результатов

nodata_value = 9999

t_xi_al = np.where(t_xi > 0, t_xi_th, t_xi)
xi_al = np.where(t_xi > 0, xi_th, xi)
xi_seasonally_thawed = np.where(t_xi <= 0, xi, 9999)
xi_seasonally_frozen = np.where(t_xi > 0, xi_th, 9999)

# Сохранение итоговых растров
write_tif(r"E:/graduate_work/v_2_results/t_xi_final_era.tif", t_xi_al, gtf, proj, nodata_value)
write_tif(r"E:/graduate_work/v_2_results/xi_final_era.tif", xi_al, gtf, proj, nodata_value)

write_tif(r"E:/graduate_work/v_2_results/xi_seasonally_frozen_final_era.tif", xi_seasonally_frozen, gtf, proj, nodata_value)
write_tif(r"E:/graduate_work/v_2_results/xi_seasonally_thawed_final_era.tif", xi_seasonally_thawed, gtf, proj, nodata_value)


