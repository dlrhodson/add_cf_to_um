#!/bin/bash
#this script extracts the files required by add_cf_to_um from code.metoffice.gov.uk
#need to use FCM and be connected to MOSRS
#Probably by being logged into PUMA!
#11/Dec/2023
#Dan Hodson

#extract the standard STASHMASTER_A file for the UM this model version
job_path0=$(grep "^job_path" cf_to_um.conf)
job_path=${job_path0##*=}
um_source0=$(grep "^config_root_path" $job_path/app/fcm_make_um/rose-app.conf)
um_source=${um_source0##*=}
echo $um_source
fcm export $um_source/rose-meta/um-atmos/HEAD/etc/stash/STASHmaster/STASHmaster_A
um_version0=$(grep UM_VERSION STASHmaster_A)
um_version=${um_version0##*=}
mv STASHmaster_A STASHmaster_A_${um_version} 
exit

#extract the common cf->stash mappings from the mip_convert tool in CDDS                                                   
fcm export fcm:cdds.x-tr/mip_convert/mip_convert/process/common_mappings.cfg
fcm export fcm:cdds.x-tr/mip_convert/mip_convert/process/HadGEM3_mappings.cfg
#extract the default UM stash from u-bg466 - GC3,1 CMIP6 Hist run                                                          
fcm export https://code.metoffice.gov.uk/svn/roses-u/b/g/4/6/6/trunk/app/um/rose-app.conf rose-app.conf.bg466
#fcm export https://code.metoffice.gov.uk/svn/um/main/branches/pkg/Share/vn10.7_CMIP6_production_mods/rose-meta/um-atmos/HEAD/etc/stash/STASHmaster/STASHmaster_A
#extract the standard NEMO defined diagnostics list                                                                        
fcm export https://code.metoffice.gov.uk/svn/nemo/branches/UKMO/dev_r5518_GO6_package/NEMOGCM/CONFIG/SHARED/field_def.xml
#field_def.xml has an error in the XML - a comment includes "<20cm" which confuses the XML parser                          
#replace "<" with "&lt;" - the correct XML representation of "<"                                                           
sed -i 's/<20cm/\&lt\;20cm/g' field_def.xml

#file for that defines all the available ice diagnostics
fcm export https://code.metoffice.gov.uk/svn/cice/main/branches/pkg/Config/vn5.1.2_GSI8.1_package_branch/cice/source/ice_history_shared.F90 
