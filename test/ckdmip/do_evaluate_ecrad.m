path(path,'../common');
figure(1)
evaluate_ckd_lw_fluxes('ckdmip_evaluation1_lw_fluxes_present_reduced.nc','ckdmip_evaluation1_lw_out.nc','ecRad','Present-day')
figure(2)
evaluate_ckd_sw_fluxes('ckdmip_evaluation1_sw_fluxes_present_reduced.nc','ckdmip_evaluation1_sw_out.nc','ecRad','Present-day')
