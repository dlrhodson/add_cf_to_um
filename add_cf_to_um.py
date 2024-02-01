#!/usr/bin/env python3

#add_cf_to_um.py
#
#Reads in a CSV file of (realm,frequency,spatial_domain)
#converts this to UM STASH, NEMO and CICE diagnostics and then adds them to an
#existing ROSE suite
#usage: add_cf_to_um.py -c <conf_file> -s <stash_type>
#conf_file is a configuration file
#stash_type is um or xios - defined which stash to use to add diagnostics for the atmosphere


from copy import deepcopy
import lxml.etree as ET
#import xml.etree.ElementTree as ET
import configparser
import argparse
import hashlib
import logging
#import uuid
import csv
import sys
import re
import os

class color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def sections_equal(section1, section2):
    # Check if both sections exist
    if config.has_section(section1) and config.has_section(section2):
        # Get key-value pairs for both sections
        items1 = config.items(section1)
        items2 = config.items(section2)

        # Check if the items in both sections are equal
        return items1 == items2
    else:
        # One or both sections do not exist
        return False


def is_equal_except(dict1,dict2,not_key):
    '''
    tests if two dicts are equal except for one key
    '''

    d1=[(key, value) for key, value in dict1 if key != not_key]
    d2=[(key, value) for key, value in dict2 if key != not_key]

    return d1==d2
    
    
def bold(message):
    return(color.BOLD+message+color.END)

def is_subset(subset_dict, main_dict):
    return all(main_dict.get(key) == value for key, value in subset_dict.items())

 
def read_cf_mappings():
    '''
    read the cf mappings from the mappings files linked to from the config file
    '''
    cf_mappings = configparser.ConfigParser()   
    # turn of the lowercaseization of keys
    cf_mappings.optionxform = lambda option: option
    #configFilePath = 'common_mappings.cfg'
    configFilePath = main_config['main']['mappings']
    if ',' in configFilePath:
       #this contains multiple config files, split into a list
       configFilePaths=configFilePath.split(',')
       for path in configFilePaths:
          if not os.path.isfile(path):
             print("ROSE conf file "+path+" does not exist")
             exit()

       cf_mappings.read(configFilePaths)
    else:
       if not os.path.isfile(configFilePath):
          print("ROSE conf file "+configFilePath+" does not exist")
          exit()
       cf_mappings.read(configFilePath)
    return(cf_mappings)


def diags_from_expression(diag):
    '''
    extracts the stash codes, nemo and cice  and cf_diags from the mapping expression for diag
    '''
    expression0=cf_mappings[diag]['expression']
    
    nemo_cice_diags=[]

    
    #    if 'umo' in diag:
    #        import pdb; pdb.set_trace()
    
    #Run through expression and extract all the UM stashcodes and other diagnostic names
   
    #remove any \n
    expression=expression0.replace('\n',' ')
    #extract any UM stash codes
    pattern_s = r'm\d{2}s\d{2}i\d{3}\[.*?\]'
    um_diags = re.findall(pattern_s, expression)
    #remove any stash codes of the form mNNsNNiNN[XXX]
    #[XXX] defines the meaning and levels etc
    expression = re.sub(pattern_s,'', expression)
    #extract all UM stash codes that do not have a following [XXX]
    pattern_s = r'm\d{2}s\d{2}i\d{3}'
    um_diags_nobracket = re.findall(pattern_s, expression)
    #remove any stash codes of the form mNNsNNiNN
    expression = re.sub(pattern_s,'', expression)
    #remove any functions of the form 'function('
    pattern=r'[0-9a-z_]+\('
    expression=re.sub(pattern,'',expression)
    #replace '*' by space - as some expressions don't include spaces around *!
    expression=expression.replace('*',' ')
    #remove commas,  + - and and ( or )
    translation_table = str.maketrans("", "", ",+-/()")
    expression=expression.translate(translation_table)
    #Now split into parts
    sub_diags=[x for x in expression.split(' ') if x !='']
    if sub_diags:
        for sub_diag in sub_diags:
            #does this string contain any lower case letters?
            lower_case=re.search(r'[a-z]',sub_diag)
            if lower_case:
                #this should contain a nemo or cice diagnostic
                if not ('=' in sub_diag or 'mask' in sub_diag or '_0' in sub_diag):
                    #if sub_diag contains an '=' this is probably a mask argument - so we can ignore
                    #if sub_diag contains 'mask' it is probably a mask - so we can ignore
                    #if sub_diag contains '_0' it is probably a REFERENCE diag from another experiment - so we will ignore it
                    #let's remove any square bracket expressions remaining (eg thetao[depth<2025])
                    sub_diag=re.sub('\[.*?\]','',sub_diag)
                    nemo_cice_diags.append(sub_diag)
    #convert to list if unique elements 
    um_diags=list(set(um_diags))
    um_diags_nobracket=list(set(um_diags_nobracket))
    nemo_cice_diags=list(set(nemo_cice_diags))
    
    return(um_diags,um_diags_nobracket,nemo_cice_diags)

def add_nemo_cice_diagnostic(diag,freq,dims):
    '''
    Adds a single diag directly to ocean or ice
    '''
        #Here the diag has no CF mapping OR we have no um diags, and only 1 nemo or cice diag and this points back to itself (e.g evs -> evs )
    #We need to look in the NEMO and CICE diagnostics for the final mappings
    if diag=='ice_present':
        #mappings in common_mappings.cfg is wrong for sitimefrac !
        #ice_present in the confing by ice_history_shared.F90 is f_icepresent
        #           fdiag='f_icepresent'
        #        else:
        #           fdiag='f_'+this_expression
        diag='icepresent'

    #is diag a cice diagnostic?
    if 'f_'+diag in cice.cice_diagnostics:
        plog(diag+" is a CICE diagnostic, adding..")
        cice.addIceDiag('f_'+diag,freq,dims)
        return()

    
    #Here diag is either a nemo diag that exists in the xml definitions or it is undefined
    #Check all the XML definitions
        
    #either diag is NOT defined in the mapping tables, OR it IS, but the definition is circular! (eg umo -> umo )
    fields=nemo.nemo_full_diagnostics.findall(".//field[@name='"+diag+"']")
    if len(fields)>0:
        #we found a matching diagnostics in the NEMO defined diagnostics
        plog("Adding NEMO diag "+diag)
        nemo.addOceanDiag(diag,freq,dims)     
        return()

    fields=nemo.nemo_diagnostic_request.findall(".//field[@name='"+diag+"']")
    if len(fields)>0:
        #we found a matching diagnostics in the NEMO user defined diagnostics
        plog("Adding NEMO diag "+diag)
        new_line=dict(line)
        new_line['variable']=diag
        nemo.addOceanDiag(diag,freq,dims)
        return()

        
    #is this diagnostics just commented out in the user defined diagnostics?
    found_in_comments=False
    for off_diag in nemo.nemo_diagnostic_request_off:
        fields=off_diag.findall(".//field[@name='"+diag+"']")
        if len(fields)==1:
            plog(diag+" found in the commented out diagnostics")
            print("uncommenting "+diag)
            
            this_freq=nemo.freq_map[freq]
            file_group=nemo.get_file_group(this_freq)
            #does this diag have a parent in the off_diag?
            parent=fields[0].getparent()
            if len(parent)==0:
               print("No parent found for "+diag+" in comments!")
               import pdb; pdb.set_trace()
            this_name_suffix=parent.attrib['name_suffix']

            #is there an EXISTING file within this filegroup that has this suffix
            existing_file=file_group.findall(".//file[@name_suffix='"+this_name_suffix+"']")
            if not existing_file:
                #no - so we need to add it
                this_id=parent.attrib['id']
                if this_id in nemo.file_element_id_list:
                    #a file with this ID already exists
                    #Need to add one
                    new_file_id=nemo.get_unique_file_id()
                    parent.attrib['id']=new_file_id
                if not 'output_freq' in parent.attrib:
                    #add a output_freq if there is not one already
                    parent.attrib['output_freq']=this_freq
                plog("Adding new file element for "+this_name_suffix)
                new_file_element=ET.SubElement(file_group,'file')
                new_file_element.attrib.update(dict(parent.attrib))

            #add this field to the file group
            print(this_name_suffix)
            #SOMETHING GOING WRONG HERE!
            nemo.add_field_to_file_group(fields[0],file_group,this_name_suffix)
            return()
            #What is the best way to uncomment this?
            #Need to go Into addOceanDiag somewhrer
        if len(fields)>1:
           print("Too many  matches of "+diag+" in the comment fields!")
           import pdb; pdb.set_trace()


    if check_output:
        nemo.nc_check_ocean(diag,freq,dims)
    #if we are here - we didn't find the diag anywhere!
    plog(diag+" not found in anywhere in NEMO diagnostics definitions")
    nemo.missing.append(diag)

                    
def add_cf_diagnostic(diag,freq,dims):
    '''
    add a cf diagnostic
    using the CF mappings (*cfg) from mipconvert
    a cf_diagnostic from a given realm (e.g. ocean) may map to a function of diagnostics from the UM, NEMO and CICE AND other cf_diagnostics
    cf_diag = function( UM_diag, NEMO_diag, CICE_diag, cf_diag_2)
    where cf_diag_2 can be ANOTHER cf_diagnostic that in turn needs mapping
    Sometimes cf mappings will map to themselves - this implies they are variables intrinsic to the model (eg NEMO or CICE)
    '''
    
    print(">>>>>>  "+diag)

    ##HERE we need to proceed if the diag is in cf_mappings
    ##IF it is, but the expression points back to itself - need to map to natice nemo or cice
    ##IF it is NOT, we need to map to native nemo or cice 
    
    
    #if not diag in cf_mappings:
    #   print("No mapping for "+diag+" in the mapping files?")
    #   import pdb; pdb.set_trace()
    #if 'evs' in diag:
    #   import pdb; pdb.set_trace()

    
    nemo_cice_diags=[]
    
    if diag in cf_mappings:
        #diag has a mapping in cf_mappings - work through each term in the mapping expression and add each sub diagnostic in that expression
        um_diags,um_diags_nobracket,nemo_cice_diags=diags_from_expression(diag)
 
        # if we have no um diags, and only 1 nemo or cice diag and this points back to itself (e.g evs -> evs )
        # then we probably have a native nemo or cice diag and need to move on to searches in the xml and conf files
        # We NOT this here for all the other cases 
        if not (len(nemo_cice_diags)==1 and nemo_cice_diags[0]==diag and len(um_diags)==0 and len(um_diags_nobracket)==0):
            #Now add the UM diags
            for um_diag in um_diags:
                if diag==um_diag:
                    #this is the should never happen!
                    import pdb; pdb.set_trace()
                #plog("Add "+um_diag)
                um.add_stash(um_diag,freq,dims)
               
                
                #WHAT ABOUT [] heres? blev and lbproc, lblev? 
                #um.add_stash(um_diag,freq,dims)

            #And the UM diags with NO post stash brackets
            for um_diag in um_diags_nobracket:
                if diag==um_diag:
                    #this is the should never happen!
                    import pdb; pdb.set_trace()
                #plog("Add "+um_diag)
                um.add_stash(um_diag,freq,dims)
               
                #um.add_stash(um_diag,freq,dims)

            #For everything else we need to recurse as these diags may have additional mappings
            for nemo_cice_diag in nemo_cice_diags:
                #is nemo_cice_diag a nemo diagnostic?
                #Skip any Ofx files
                if nemo_cice_diag in cf_mappings:
                    if cf_mappings[nemo_cice_diag]['mip_table_id']=='Ofx':
                        plog(nemo_cice_diag+" is just an Ofx field- skipping")
                        continue
                        #skip to next nemo_cice_diag
                if not nemo_cice_diag==diag:
                   #only recurse IF we are not about to enter an infinite loop!
                   # eg  evs -> func (a,b,evs)
                   
                   add_cf_diagnostic(nemo_cice_diag,freq,dims)
                else:
                   plog("Nemo diag and diag matches!")
                   plog("Adding "+diag+" directly")
                   add_nemo_cice_diagnostic(diag,freq,dims)
                   

            return()

    #Here the diag has no CF mapping OR we have no um diags, and only 1 nemo or cice diag and this points back to itself (e.g evs -> evs )
    #We need to look in the NEMO and CICE diagnostics for the final mappings

    if nemo_cice_diags:
        #we have a single nemo or cice diagnostic that exists in the cf mappings
        if len(nemo_cice_diags)>1:
            print("Something went wrong!")
            import pdb; pdb.set_trace()

   
    #Here nemo_cice_diag is either [] OR nemo_cice_diag == diag
    #so, use diag from here on

    add_nemo_cice_diagnostic(diag,freq,dims)
    return()
            
            

 
#############  Atmosphere/Land class for UM
class UM:
    #UM stash  class
    # add_cf_diagnostic() adds an atmosphere/land/landice cf variable as the require STASH codes
    def __init__(self,umOrXIOS):
        #we can use the STASH in the um app or the STASH in the xml app (for netcdf)
        #umOrXIOS = 'um' or 'xios'
        #reads in all configuration files
        #and sets up all mappings

        self.nc_found=[] #list of UM diagnostics found in the NC file during --check_output
        self.nc_missing=[] #list of UM diagnostics missing in the NC file during --check_output
        self.missing=[] # list of diagnostics we failed to add!
        self.added=[] # list of diagnostics we succeeded in adding!
        self.umOrXIOS=umOrXIOS
        self.default_usage={'mon':"'UPM'",
                            'day':"'UPD'"
                            }

        #Mappings between LBPROC and ityp in time domains
        self.lbproc_mappings={'128':'3',
                              '4096':'6',
                              '8192':'5'
                              }
        #defining day and month in the TIME DOM
        #day is ever 1 days (3) sampling every time step (1)
        #mon is every 30 days (3), sampling every time step (1)
        self.tim_dom={'day':{'ifre': '1', 'unt1': '3', 'unt2': '1', 'unt3': '3'},
                      'mon':{'ifre': '30', 'unt1': '3', 'unt2': '1', 'unt3': '3'},
                      '6hrPt':{'ifre': '6', 'istr':'6', 'ityp':'1', 'unt3': '2'}
                      }
        
        
        self.freq_mappings={'mon':"'TMONMN'",
                            'day':"'TDAYMN'"
                            }
        self.space_mappings={'longitude latitude time':"'DIAG'",
                             'longitude latitude height2m time':"'DIAG'",
                             'longitude latitude height10m time':"'DIAG'",
                             'longitude latitude typesi time':"'DIAG'",  
                             'longitude latitude alevel time':"'DALLTH'",
                             'longitude latitude alevhalf time':"'DALLRH'",
                             'longitude latitude plev8 time':"'PLEV8'",
                             'longitude latitude plev19 time':"'PLEV19'",
                             'longitude latitude plev14 time':"'PLEV14'",
                             'longitude latitude soil time':"'DSOIL'",
                             'longitude latitude sdepth time':"'DSOIL'",
                             'longitude latitude sdepth1 time':"'DSOIL1'"
                             }

        self.xios_stream_ids=[]
        self.xios_stream_filename_bases=[]
        self.cmip6_use_mappings={}
        self.use_list=[]
        self.cf_to_stash={}
        self.rose_time_domain_mappings={}
        self.rose_space_domain_mappings={}
        self.use_matrix={}
        self.stash_levels={}
        self.stash_names={}
        self.stash_pseudo_levels={}  #pseudo level mapping for stash codes
        #Stash mapping
        self.stashMap = configparser.ConfigParser()  
        #turn of the lowercaseization of keys
        self.stashMap.optionxform = lambda option: option
        self.configFilePath = main_config['main']['mappings']
        #stashMap.read(configFilePath)
        self.read_stash_mappings()
        self.read_STASHmaster_A_levels()
        #main rose-app.conf for this Job
        valid_options=['um','xios']
        if umOrXIOS in valid_options:
            if umOrXIOS=='xios':
                self.rose_stash=main_config['user']['job_path']+'app/xml/rose-app.conf'
            else:
                self.rose_stash=main_config['user']['job_path']+'app/um/rose-app.conf'
        else:
            print(umOrXIOS+" is not a valid UM class type")
            import pdb; pdb.set_trace()

        self.rose,self.rose_header=self.read_rose_app_conf(self.rose_stash)
        #rose,rose_header=self.read_rose_app_conf(rose_stash)
        #CMIP6 map


        if umOrXIOS=='xios': 
           #read in all the existing xios stream definitions
           self.get_xios_streams()
        else:
           print("UM- need to add the USAGE stream handler!")
           import pdb; pdb.set_trace()


        cmip6_rose=main_config['main']['cmip6']
        #read in standard CMIP6 time and spatial domain definitions
        self.cmip6,self.cmip6_header=self.read_rose_app_conf(cmip6_rose)

        #get list of usages in ROSE -> use_list
        self.get_use_list()
        #get a dict of mappings of UPX->umstash_use()
        self.get_cmip6_use_mappings()
        
        #matrix mapping time and space to USAGE for output
        self.get_use_matrix()

        #create mappings for time domain cmip6 -> rose
        #use the initial mappings for cf_time -> stash time domain label in freq_mappings
        #find this domain label in the cmip6 reference, then check to see if any domain usages in rose is IDENTICAL (in terms of namelist)
        #If so, use THAT domain label (dom_name) from ROSE
        #otherwise we need to ADD the cmip6 domain to ROSE
        #self.rose_time_domain_mappings=get_time_mappings(freq_mappings)
        self.get_time_mappings()

        #create mappings for space domain cmip6 -> rose
        #as for time domain, but for space domin
        #self.rose_space_domain_mappings=get_space_mappings(space_mappings)
        self.get_space_mappings()



    def copy_cmip6_tim_dom_to_rose(self,this_cmip6,this_domain):
        '''
        copy a cmip6 time domain to the rose stash
        '''

        cmip6_time=this_cmip6.name
        cmip6_tim_dom={}
        cmip6_tim_dom_full={}
        for key in this_cmip6:
           this_item=this_cmip6[key]
           cmip6_tim_dom_full[key]=this_item
           if not (('!!' in key) or ('meta' in key)):
              cmip6_tim_dom[key]=this_item
           
        plog("Adding "+this_cmip6['tim_name']+" across from the CMIP6 reference to the "+bold(self.umOrXIOS)+" STASH")
        #copy across CMIP6 TIME DOMAIN
        if self.umOrXIOS=='xios':
            #this is the XIOS STASH - we need to add the ts_enabled=.false. item
            cmip6_tim_dom_full['ts_enabled']='.false.'

        #the uuid hash for this time may not be the same as if it had been calculated for this version of the UM
        #so let's recalculate it:
        this_uuid=self.get_uuid_hash(cmip6_tim_dom_full)
        cmip6_time_sub=cmip6_time.split('(')
        cmip6_time_new=cmip6_time_sub[0]+'('+cmip6_time_sub[1].split('_')[0]+'_'+this_uuid+')'
        self.rose[cmip6_time_new]=cmip6_tim_dom_full
        #import pdb; pdb.set_trace()

        self.rose_time_domain_mappings[freq+'_'+this_cmip6['ityp']]=this_cmip6['tim_name']
        #Now need to add this Time domain to the use_matrix
        plog("Now need to add "+this_cmip6['tim_name']+" to the use matrix")
        #find all variables in cmip6 that use this time domain
        cmip6_variable_keys=[key for key in self.cmip6.keys() if 'umstash_streq' in key]
        #does any existing cmip6 stash diag using this tim_name?
        this_time=this_cmip6['tim_name']
        cmip6_time_match=[key for key in cmip6_variable_keys if self.cmip6[key]['tim_name']==this_time]
        if not cmip6_time_match:
            plog("No existing cmip6 diag uses "+this_cmip6['tim_name']+" !")
            this_use=self.default_usage[freq]
            plog("So we will use the default usage for "+freq+" : "+this_use)
            #What domain?
        else:
            this_use=self.cmip6[cmip6_time_match[0]]['use_name']
            this_domain=self.cmip6[cmip6_time_match[0]]['dom_name']
            print("Found existing use of "+this_time+" which uses "+this_use+" and "+this_domain)
            #import pdb; pdb.set_trace()

            

        if not this_time in self.use_matrix:
           #time not already in use matrix - add
           self.use_matrix[this_time]={this_domain:[this_use]}
           self.check_use_already_exists(this_use)
                 #if not this_use in self.use_list:
                 #this usage does not exist in ROSE
                 #   print("FLOSH")
                 #   print("copying usage "+this_use+" to ROSE")
                 #   #copy usage with name this_use from CMIP6
                 #   new_use=self.cmip6[self.cmip6_use_mappings[this_use]]
                 #   #create a consistent new file id for this usage
                 #   new_use['file_id']=self.create_new_file_id()
                 #   #copy reference to this usage to ROSE
                 #   self.rose[self.cmip6_use_mappings[this_use]]=new_use
                 #   
                 #   import pdb; pdb.set_trace()
                 #   #####HEEEREE
        else:#
            #this_time is already in the use_matrix
            if this_domain in self.use_matrix[this_time]:
                #this domain already exist here
                if not this_use in self.use_matrix[this_time][this_domain]:
                    #only add if this use isn't already in the list
                    self.use_matrix[this_time][this_domain].append(this_use)
                    check_use_already_exists(this_use)
                    if not this_use in self.use_list:
                        print("FISH")
                        #this usage does not exist in ROSE
                        print("copying usage "+this_use+" to ROSE")
                        import pdb; pdb.set_trace()
                        #####HEEEREE
            else:
                #domain not here yet
                self.use_matrix[this_time][this_domain]=[this_use]
                #check to ensure that this_use exists and points to an output stream - if not TAKE ACTION!
                self.check_use_already_exists(this_use)
                if not this_use in self.use_list:
                    print("HELLO")
                    #this usage does not exist in ROSE
                    plog("copying usage "+this_use+" to ROSE")
                    #copy usage with name this_use from CMIP6
                    new_use=self.cmip6[self.cmip6_use_mappings[this_use]]
                    #create a consistent new file id for this usage
                    new_use['file_id']=create_new_file_id()
                    #copy reference to this usage to ROSE
                    self.rose[self.cmip6_use_mappings[this_use]]=new_use
                    import pdb; pdb.set_trace()
                    #####HEEEREE
        plog(this_time+" added to Use Matrix")




        
    def custom_sort(self,item):
       #custom sort for secton items that should ignore the prefix "!!"
       return(item.lstrip('!!'))

        #Now rose should have all the required time and space domains defined as in the freq_mappings and space_mappings
    def get_uuid_hash(self,section):
       #computes the correct hash for this stash (as in TidyStashValidate in stash_indices.py)
       #NOTE - this currently ONLY works for streq - use and domain will need to have the name fields removed before this is done!
       text=''
       
       if 'isec' in section:
          #this must be a streq:
          section_keys=[key for key in section]
       else:
          #must be domain or use or something?
          #if so, we want to exclude the NAME label from the hash calclation
          #this allows us to compare sections that are identical, except for the names
          section_keys=[key for key in section if not key in ['use_name','dom_name','tim_name']]
       #sort the keys, ignoring the '!!' prefix during the sort
       section_keys.sort(key=self.custom_sort)
       for key in section_keys:
          #sometimes the value can come from a multi-line and may have \n and = - remove these for the uuid hash
          this_value=str(section[key]).replace('\n','').replace('=','')
          text+=key+'='+this_value+'\n'

       #if 'rlev' in text:
       #   import pdb; pdb.set_trace()

       uuid=hashlib.sha1(text.encode(encoding="utf8")).hexdigest()[:8]
     
       return(uuid)


       
    def check_use_already_exists(self,use):
       
       #check to ensure that this_use exists and points to an output stream - if not TAKE ACTION!
       if not use in self.use_list:
          plog(use+" does not already exist in ROSE")
          if not 'usage' in main_config:
             print("No [usage] section in "+conf_file)
             print("Please add the following to "+conf_file+" and then re-run ")
             print()
             print("[usage]")
             print(use.replace("'","")+"=<xios_stream_reference>")
             print()
             print("Where <xios_stream_reference> is either one of the following existing xios streams:")
             for i in self.xios_stream_ids:
                      print(i)
             print()
             print("Or a new xios_stream that you want to define")
             exit()
          else:
             usage_section=main_config['usage']
             #remove ' just in case
             if not use.replace("'","") in usage_section:
                print("[usage] section exists in "+conf_file+ " but no usage stream mapping defined for "+use+" in this [usage] section")
                print("Please add a use mapping to the [usage] section - something like:")
                print(use.replace("'","")+"=<xios_stream_reference>")
                print()
                print("Where <xios_stream_reference> is either one of the following existing xios streams:")
                for i in self.xios_stream_ids:
                   print(i)
                print()
                print("Or a new xios_stream that you want to define")
                exit()
             else:
                
                #OK so we have a usage mapping!
                #remove ' just in case
                use_stream=usage_section[use.replace("'","")]
                #does this stream already exist?
                if not use_stream in self.xios_stream_ids:
                   plog(use_stream+" does not exist in the existing output streams ")

                   #do we have this stream defined in the confi?
                   xios_config=[key for key in main_config if 'xios_streams' in key]
                   xios_config_map={}
                   for key in xios_config:
                      xios_config_map[main_config[key]['file_id'].replace("'","")]=key
                   if use_stream in xios_config_map:
                      plog("Aha - I see we have "+use_stream+" defined in "+conf_file)
                      #check we don't already use the filename_base
                      this_xios_stream=xios_config_map[use_stream]
                      this_filename_base=main_config[this_xios_stream]['filename_base'].replace("'","")

                      if this_filename_base in self.xios_stream_filename_bases:
                         print("Oh - the filename base name for the new xios_stream ("+this_filename_base+") is already used by another xios stream!")
                         print("Existing stream base names are:")
                         for i in self.xios_stream_filename_bases:
                            print(i)
                         print("Please change this and rerun")
                         exit()
                      else:
                         #Everthing OK - we now need to add the XIOS stream to ROSE
                         plog("Adding "+use_stream+" to ROSE")
                         self.add_new_xios_stream(this_xios_stream)
                         #and the USAGE
                         plog("Adding "+use+" to ROSE")
                         self.add_new_usage(use,use_stream)
                         print("DONE")
                         return()
                         
                   else:
                      for i in self.xios_stream_ids:
                         print(i)
                      
                      print("Please either adjust this to an existing output stream - or add an xios stream definition section to "+conf_file+". Something like")
                      print()
                      print("[namelist:xios_streams("+use_stream+")]")
                      print("compression_level = 0")
                      print("file_id = '"+use_stream+"'")
                      print("""!!filename = ''
filename_base = './${RUNID}a_mon_'
l_reinit = .true.
output_freq_unit = 4
output_freq_value = 1
reinit_end = -1
reinit_start = 0
reinit_step = 1
reinit_unit = 4
time_counter = centered_exclusive
time_stamp_name = 'creation_date'
timeseries = none
!!ts_prefix = ''
!!use_ts_prefix = .false.
uuid_name = 'tracking_id'
                   """)
                      exit()
                else:
                   plog("using mapping for "+use+" in "+conf_file)
                   plog(use+"->"+use_stream)
                   #use stream mapping exists! Add New usage and stream
                   self.add_new_usage(use,use_stream)
                   return()
                import pdb; pdb.set_trace()
       else:
          #this use already exists in the use list - do nothing
          return()
       
    def add_new_xios_stream(self,xios):
       #link the new xios stream into ROSE
       self.rose[xios]=main_config[xios]


          
    def add_new_usage(self,use,use_stream):
       #creates a new use mapping and adds to ROSE
       plog("Creating new usage "+use+" using "+use_stream)
       if not "'" in use_stream:
          #use_stream should be e.g "'xios_upm_1m'"
          #add single quote around, if required
          use_stream="'"+use_stream+"'"
       new_use={'file_id':use_stream,'locn':3,'macrotag':0,'use_name':use}
       #compute the uuid hash for this usage
       hex_uuid=self.get_uuid_hash(new_use)
       umstash_use="namelist:umstash_use("+use.replace("'","").lower()+"_"+hex_uuid+")"
       self.rose[umstash_use]=new_use
       

       #hex_uuid = format(int(uuid.uuid4().hex[:8], 16), 'x')
       
       plog("Adding to use_list")
       self.use_list.append(use)
       print("Done")

          
    def get_xios_streams(self):
       #extracts all the xios_stream file_id s to self.xios_strea_ids
       xios_stream_keys=[key for key in self.rose.keys() if 'xios_streams' in key]
       for key in xios_stream_keys:
          xios_stream=self.rose[key]
          xios_stream_name=xios_stream['file_id'].replace("'","")
          xios_stream_filename_base=xios_stream['filename_base'].replace("'","")
          self.xios_stream_ids.append(xios_stream_name)
          self.xios_stream_filename_bases.append(xios_stream_filename_base)
          
       
       
           
    def create_new_file_id(self):
       #creates a new File ID for STASH output
       #read in existing files output file stream and assign a new one
       #or many read in definition from the config file?
       #####HERE
       #where are the current list fo XIOS files held?
       #xml/rose-app.conf - I think
       #xios_streams
       import pdb; pdb.set_trace()

       return("'NEW'")
    

    def read_STASHmaster_A_levels(self):
        file=main_config['main']['stashmaster_A']
        if not os.path.isfile(file):
            print(file+" does not exist")
            exit()

        stashfile=open(file,'r')
        stashm=stashfile.readlines()
        #https://reference.metoffice.gov.uk/um/c4/_level_type_code
        #https://code.metoffice.gov.uk/doc/um/vn13.4/papers/umdp_C04.pdf    
        #1. Model rho levels
        #2. Model theta levels
        #3. Pressure levels
        #4. Geometric height levels
        #5. Single level
        #6. Deep soil levels
        #7. Potential temperature levels
        #8. Potential vorticity levels
        #9. Cloud threshold levels
        level_names={0:'unspecified',
                     1:'atm_rho',
                     2:'atm_theta',
                     3:'pressure',
                     4:'unknown',
                     5:'single',
                     6:'soil',
                     7:'theta'
                     }

        for line in stashm:
            if '|' in line:
               #see https://code.metoffice.gov.uk/doc/um/vn13.4/papers/umdp_C04.pdf
                bits=line.split('|')
                #line 1 in stash master
                #|Model |Sectn | Item |Name 
                if bits[0]=='1':
                    model=bits[1].strip(' ')
                    if model=='-1':
                        #end of file
                        break
                    sec=bits[2].strip(' ')
                    item=bits[3].strip(' ')
                    name=bits[4]
                    scode='m0'+model+'s'+sec.zfill(2)+'i'+item.zfill(3)
                    self.stash_names[scode]=name
                #line 2
                #|Space |Point | Time | Grid |LevelT|LevelF|LevelL|PseudT|PseudF|PseudL|LevCom|
                #  1       2       3      4       5    6      7      8     9       10     11
                if bits[0]=='2':
                    headings=['Space','Point','Time','Grid','LevelT','LevelF','LevelL','PseudT','PseudF','PseudL','LevCom']
                    this_line={}

                    
                    for bit in range(len(headings)):
                       this_line[headings[bit]]=int(bits[bit+1].strip(' '))
                    self.stash_levels[scode]=this_line
                    if not this_line['LevelT'] in level_names:
                        print("Unknown level! "+str(this_line['LevelT']))
                        print(name)
                        import pdb; pdb.set_trace()
                    else:
                        self.stash_levels[scode]=this_line
        return()



     
    def read_STASHmaster_A_levels_old(self):
        file=main_config['main']['stashmaster_A']
        if not os.path.isfile(file):
            print(file+" does not exist")
            exit()

        stashfile=open(file,'r')
        stashm=stashfile.readlines()
        #https://reference.metoffice.gov.uk/um/c4/_level_type_code
        #https://code.metoffice.gov.uk/doc/um/vn13.4/papers/umdp_C04.pdf    
        #1. Model rho levels
        #2. Model theta levels
        #3. Pressure levels
        #4. Geometric height levels
        #5. Single level
        #6. Deep soil levels
        #7. Potential temperature levels
        #8. Potential vorticity levels
        #9. Cloud threshold levels
        level_names={0:'unspecified',
                     1:'atm_rho',
                     2:'atm_theta',
                     3:'pressure',
                     4:'unknown',
                     5:'single',
                     6:'soil',
                     7:'theta'
                     }

        for line in stashm:
            if '|' in line:
               #see https://code.metoffice.gov.uk/doc/um/vn13.4/papers/umdp_C04.pdf
                bits=line.split('|')
                #line 1 in stash master
                #|Model |Sectn | Item |Name 
                if bits[0]=='1':
                    model=bits[1].strip(' ')
                    if model=='-1':
                        #end of file
                        break
                    sec=bits[2].strip(' ')
                    item=bits[3].strip(' ')
                    name=bits[4]
                    scode='m0'+model+'s'+sec.zfill(2)+'i'+item.zfill(3)
                #line 2
                #|Space |Point | Time | Grid |LevelT|LevelF|LevelL|PseudT|PseudF|PseudL|LevCom|
                #  1       2       3      4       5    6      7      8     9       10     11
                if bits[0]=='2':
                    level=int(bits[5].strip(' '))
                    if not level in level_names:
                        print("Unknown level! "+str(level))
                        print(name)
                        #import pdb; pdb.set_trace()
                    else:
                        self.stash_levels[scode]=level
                    pseudo_level=int(bits[8].strip(' '))
                    #pseudo levels, eg SW rad bands, surface tile types, sea ice categories
                    self.stash_pseudo_levels[scode]=pseudo_level
        return()


    def read_stash_mappings(self):
        self.cf_to_stash = configparser.ConfigParser()   
        #turn of the lowercaseization of keys
        self.cf_to_stash.optionxform = lambda option: option
        #configFilePath = 'common_mappings.cfg'
        configFilePath = main_config['main']['mappings']
        if ',' in configFilePath:
           #this contains multiple config files, split into a list
           configFilePaths=configFilePath.split(',')
           for path in configFilePaths:
              if not os.path.isfile(path):
                 print("ROSE conf file "+path+" does not exist")
                 exit()

           self.cf_to_stash.read(configFilePaths)
        else:
           if not os.path.isfile(configFilePath):
              print("ROSE conf file "+configFilePath+" does not exist")
              exit()
           self.cf_to_stash.read(configFilePath)
        return()


        
    def read_rose_app_conf(self,file):
        config = configparser.ConfigParser() 
        config.optionxform = lambda option: option
        
        configFilePath = file
        #sometimes the first line is NOT a section header, but some other string!
        #test if this is true, if so, extract this header
        header_flag=True
        if not os.path.isfile(file):
            print("ROSE conf file "+file+" does not exist")
            exit()
 
        with open(configFilePath) as stream:
            header,config_string=stream.read().split('\n',1)
            
            if '[' in header:
                 config.read_string(header+'\n'+config_string)
                 header=''
            else:
                config.read_string(config_string)
            #config.read_string("[DUMMY]\n" + stream.read())
        return(config,header)


    def get_use_list(self):
        rose_use_keys=[key for key in self.rose.keys() if 'umstash_use' in key]
        for key in rose_use_keys:
            self.use_list.append(self.rose[key]['use_name'])

    def get_cmip6_use_mappings(self):
        cmip6_use_keys=[key for key in self.cmip6.keys() if 'umstash_use' in key]
        for key in cmip6_use_keys:
            self.cmip6_use_mappings[self.cmip6[key]['use_name']]=key
                                    
        
    def get_time_mappings(self):
        #loop over all freq_mappings
        #get all time domain keys for cmip6
        #create a mapping for the rose time domain
        #WILL NEED TO ADD ANY MISSING TIME DOMAINS TO ROSE
        cmip6_time_keys=[key for key in self.cmip6.keys() if 'umstash_time' in key]
        rose_time_keys=[key for key in self.rose.keys() if 'umstash_time' in key]

        #rose_freq_mappings={}
        for freq in self.freq_mappings:
            this_time_dom=self.freq_mappings[freq]
            #loop over all cmip6 time domains
            tim_name_found=False
            for cmip6_time in cmip6_time_keys:
                this_cmip6=self.cmip6[cmip6_time]
                #find time domain with same tim_name in cmip6
                if this_time_dom==this_cmip6['tim_name']: #.replace("'",""):
                    tim_name_found=True
                    cmip6_tim_dom={}
                    cmip6_tim_dom_full={}
                    for key in this_cmip6:
                        this_item=this_cmip6[key]
                        cmip6_tim_dom_full[key]=this_item
                        if not (('!!' in key) or ('meta' in key)):
                            cmip6_tim_dom[key]=this_item
                    #time domain found in cmip6 - leave loop
                    break
            if not tim_name_found:
                print('Time domain '+this_time_dom+" Not found in CMIP6 reference!")
                print("Not sure what to do now!")
                import pdb; pdb.set_trace()

            #does one of the time domains in rose match cmip6_tim_dom?
            rose_tim_found=None
            for key in rose_time_keys:
                this_rose=self.rose[key]
                rose_tim_dom={}
                #build dict of active key value pairs
                for item in this_rose:
                    this_item=this_rose[item]
                    #XML Time Domain definitions have an extra field: ts_enabled - this is from XIOS
                    #https://forge.ipsl.jussieu.fr/ioserver/raw-attachment/wiki/WikiStart/XIOS_reference_guide.pdf
                    #usually set to .false.
                    #we will ignore this for the comparisons - but need to insert if we copy accross from CMIP6 ref
                    if not (('!!' in item) or ('meta' in item) or ('ts_enabled' in item)):
                        rose_tim_dom[item]=this_item
                if rose_tim_dom==cmip6_tim_dom:
                    rose_tim_found=key
                    break
            if rose_tim_found==None:
                print("Time domain for "+freq+" not found in Rose")
                print("Need to copy across from the CMIP6 reference")
                plog("Adding "+this_cmip6['tim_name']+" across from the CMIP6 reference to the "+bold(self.umOrXIOS)+" STASH")
                #copy across CMIP6 TIME DOMAIN
                if self.umOrXIOS=='xios':
                    #this is the XIOS STASH - we need to add the ts_enabled=.false. item
                    cmip6_tim_dom_full['ts_enabled']='.false.'

                #the uuid hash for this time may not be the same as if it had been calculated for this version of the UM
                #so let's recalculate it:
                this_uuid=self.get_uuid_hash(cmip6_tim_dom_full)
                cmip6_time_sub=cmip6_time.split('(')
                cmip6_time_new=cmip6_time_sub[0]+'('+cmip6_time_sub[1].split('_')[0]+'_'+this_uuid+')'
                self.rose[cmip6_time_new]=cmip6_tim_dom_full
                #import pdb; pdb.set_trace()
                #time domain mapping is set to the NAME (eg TDAYMN) PLUS the ITYP - eg 3:mean 5:min etc
                self.rose_time_domain_mappings[freq+'_'+this_cmip6['ityp']]=this_cmip6['tim_name']
                #Now need to add this Time domain to the use_matrix
                print("Now need to add "+this_cmip6['tim_name']+" to the use matrix")
                #find all variables in cmip6 that use this time domain
                cmip6_variable_keys=[key for key in self.cmip6.keys() if 'umstash_streq' in key]
                for cmip_variable_key in cmip6_variable_keys:
                   this_cmip6_variable=self.cmip6[cmip_variable_key]
                   this_time=this_cmip6_variable['tim_name']
                   if this_time==this_cmip6['tim_name']:
                      this_domain=this_cmip6_variable['dom_name']
                      this_use=this_cmip6_variable['use_name']
                   
                      if not this_time in self.use_matrix:
                         #time not already in use matrix - add
                         self.use_matrix[this_time]={this_domain:[this_use]}
                         self.check_use_already_exists(this_use)
                         #if not this_use in self.use_list:
                         #this usage does not exist in ROSE
                         #   print("FLOSH")
                         #   print("copying usage "+this_use+" to ROSE")
                         #   #copy usage with name this_use from CMIP6
                         #   new_use=self.cmip6[self.cmip6_use_mappings[this_use]]
                         #   #create a consistent new file id for this usage
                         #   new_use['file_id']=self.create_new_file_id()
                         #   #copy reference to this usage to ROSE
                         #   self.rose[self.cmip6_use_mappings[this_use]]=new_use
                         #   
                         #   import pdb; pdb.set_trace()
                         #   #####HEEEREE
                      else:#
                         if this_domain in self.use_matrix[this_time]:
                            #this domain already exist here
                            if not this_use in self.use_matrix[this_time][this_domain]:
                               #only add if this use isn't already in the list
                               self.use_matrix[this_time][this_domain].append(this_use)
                               check_use_already_exists(this_use)
                               if not this_use in self.use_list:
                                  print("FISH")
                                  #this usage does not exist in ROSE
                                  print("copying usage "+this_use+" to ROSE")
                                  import pdb; pdb.set_trace()
                                  #####HEEEREE
                         else:
                            #domain not here yet
                            self.use_matrix[this_time][this_domain]=[this_use]
                            #check to ensure that this_use exists and points to an output stream - if not TAKE ACTION!
                            self.check_use_already_exists(this_use)
                            if not this_use in self.use_list:
                               print("HELLO")
                               #this usage does not exist in ROSE
                               print("copying usage "+this_use+" to ROSE")
                               #copy usage with name this_use from CMIP6
                               new_use=self.cmip6[self.cmip6_use_mappings[this_use]]
                               #create a consistent new file id for this usage
                               new_use['file_id']=create_new_file_id()
                               #copy reference to this usage to ROSE
                               self.rose[self.cmip6_use_mappings[this_use]]=new_use
                               import pdb; pdb.set_trace()
                               #####HEEEREE
                plog(this_time+" added to Use Matrix")


            else:
                #print("Time domain found in Rose")
                plog(self.rose[rose_tim_found]['tim_name']+" found in Rose")

                self.rose_time_domain_mappings[freq+'_'+self.rose[rose_tim_found]['ityp']]=self.rose[rose_tim_found]['tim_name']
        #return(rose_freq_mappings)
        return()


    
    
    def get_space_mappings(self):
        #loop over all space_mappings
        #get all space domain keys for cmip6
        #create a mapping for the rose space domain
        #WILL NEED TO ADD ANY MISSING SPACE DOMAINS TO ROSE
        cmip6_space_keys=[key for key in self.cmip6.keys() if 'umstash_domain' in key]
        rose_space_keys=[key for key in self.rose.keys() if 'umstash_domain' in key]
 
        #rose_space_mappings={}
        #loop over all the space mappings (e.g. 'longitude latitude time':"'DIAG'")
        #make sure that these spatial domains exist in rose, if not -copy across
        for space in self.space_mappings:
            
            this_space_dom=self.space_mappings[space]
            #loop over all cmip6 space domains
            space_name_found=False

            #loop over all the spatial domains in CMIP6
            for cmip6_space in cmip6_space_keys:
                this_cmip6=self.cmip6[cmip6_space]
                #print(this_cmip6['dom_name'])

                #find space domain with same space_name in cmip6
                #remove ' from inside cmip6 string
                if this_space_dom==this_cmip6['dom_name']: #.replace("'",""):
                    space_name_found=True
                    cmip6_space_dom={}
                    cmip6_space_dom_full={}
                    for key in this_cmip6:
                        this_item=this_cmip6[key]
                        cmip6_space_dom_full[key]=this_item
                        if not (('!!' in key) or ('meta' in key)):
                            cmip6_space_dom[key]=this_item
                    break
            if not space_name_found:
                print('Space domain '+this_space_dom+" Not found in CMIP6 reference!")
                user_domains=[key for key in main_config if 'umstash_domain' in key]
                if user_domains:
                    print("Found some user defined domains "+' '.join(user_domains))
                    user_domain_found=[key for key in user_domains if main_config[key]['dom_name']==this_space_dom]
                    if not user_domain_found:
                        print("No "+this_space_dom+" found in the user domaines :(")
                        import pdb; pdb.set_trace()
                    #just pick first - they all match the this_space_dom!
                    user_space=user_domain_found[0]
                    user_space_dom_full=main_config[user_space]
                    #copy across USER DOMAON
                    #the uuid hash for this domain may not be the same as if it had been calculated for this version of the UM
                    #so let's recalculate it:
                    this_uuid=self.get_uuid_hash(user_space_dom_full)
                    user_space_sub=user_space.split('(')
                    user_space_new=user_space_sub[0]+'('+user_space_sub[1].split('_')[0]+'_'+this_uuid+')'
                    self.rose[user_space_new]=user_space_dom_full
                    self.rose_space_domain_mappings[space]=user_space_dom_full['dom_name']
                    print("Copied user domain "+user_space_dom_full['dom_name']+" to ROSE")
                    
                    continue

                print("Not sure what to do now!")
                import pdb; pdb.set_trace()

            #does one of the space domains in rose match cmip6_tim_dom?
            rose_space_found=None
            for key in rose_space_keys:
                this_rose=self.rose[key]
                rose_space_dom={}
                #build dict of active key value pairs
                for item in this_rose:
                    this_item=this_rose[item]
                    if not (('!!' in item) or ('meta' in item)):
                        rose_space_dom[item]=this_item
                if rose_space_dom==cmip6_space_dom:
                    rose_space_found=key
                    break
            if rose_space_found==None:
                print("Space domain "+this_cmip6['dom_name']+"not found in Rose")
                plog("copying "+this_cmip6['dom_name']+" across from the CMIP6 reference")
                #logging.info("Adding "+this_cmip6['dom_name']+" across from the CMIP6 reference")
                #copy across CMIP6 DOMAON
                #the uuid hash for this domain may not be the same as if it had been calculated for this version of the UM
                #so let's recalculate it:
                this_uuid=self.get_uuid_hash(cmip6_space_dom_full)
                cmip6_space_sub=cmip6_space.split('(')
                cmip6_space_new=cmip6_space_sub[0]+'('+cmip6_space_sub[1].split('_')[0]+'_'+this_uuid+')'
                self.rose[cmip6_space_new]=cmip6_space_dom_full
                self.rose_space_domain_mappings[space]=this_cmip6['dom_name']
            else:
                #print("Space domain found in Rose")
                #logging.info(self.rose[rose_space_found]['dom_name']+" found in ROSE")
                plog(self.rose[rose_space_found]['dom_name']+" found in ROSE")
                self.rose_space_domain_mappings[space]=self.rose[rose_space_found]['dom_name']
        return()


    def get_use_matrix(self):
        #this returns a dict (matrix) use_matrix[time][space]=[usage_list]
        #usage_list is a list of usage labels that are used for that time and space pairs in ROSE use
        stash_c= [key for key in self.rose.sections() if re.search('umstash_streq', key)]
        use_mappings_t={}
        use_mappings_s={}
        #use_matrix={}

        for key in stash_c:
            this_stash=self.rose[key]
            this_use=this_stash['use_name']
            this_time=this_stash['tim_name']
            this_space=this_stash['dom_name']
            tmp_dom={}
            tmp_dom[this_space]=[this_use]

            if not this_time in self.use_matrix:
                #this_time has not yet been added to use_matrix
                self.use_matrix[this_time]=tmp_dom
            else:
                if not this_space in self.use_matrix[this_time]:
                    #this_time exists in use_matrix, but there is no entry for this_space
                    self.use_matrix[this_time][this_space]=[this_use]
                else:
                    #this_time and this_space exist in use_matrix already - is this a repeated entry for this_use, or an additional use?
                    if not this_use in self.use_matrix[this_time][this_space]:
                        #this is an additional use for this element:
                        self.use_matrix[this_time][this_space].append(this_use)
        return()


    def get_stash_codes(self,cf_variable):
        

        if not cf_variable in self.cf_to_stash.sections():
            print(cf_variable+" not found in the stash mapping!")
            plog(bold("skipping "+cf_variable))
            self.missing.append(cf_variable)
            return(None,None)
                 
        this_map=self.cf_to_stash[cf_variable]
        this_expression=this_map['expression']
        this_dimension=this_map['dimension']
        pattern = r'm\d{2}s\d{2}i\d{3}\[.*?\]'
        matches = re.findall(pattern, this_expression)
        stash_codes={}
        for match in matches:
            stash_code,attributes=match.replace(']','').split('[')
            stash_attributes={}
            for attribute in attributes.split(','):
                key,value=attribute.split('=')
                stash_attributes[key]=value
            stash_codes[stash_code]=stash_attributes
       
        return(stash_codes,this_dimension)


    def add_cf_diagnostic(self,line):
        #import pdb; pdb.set_trace()
        cf_var=line['variable']
        freq=line['time']
        dims=line['space']

        print(cf_var+" "+freq)

        #get the stash code and dimensions for this cf variable
        req_stash,dimensions=self.get_stash_codes(cf_var)

        if req_stash==None:
            #that didn't work!
            return()
        if dimensions!=dims:
            print("dimensions clash!")
            #plog("CF mappings requests ["+dimensions+"] but input diagnostics request wants ["+dims+"]")
            #plog("Using ["+dims+"]")
            plog("mappings requests ["+dimensions+"] but input diagnostics request wants ["+dims+"]\nUsing ["+dims+"]")
            plog("Using ["+dims+"]")
            dimensions=dims

        print(req_stash)
        import pdb; pdb.set_trace()

        if not freq in self.rose_time_domain_mappings:
            print("\t"+freq+" does not exist in time mappings?")
            import pdb; pdb.set_trace()

        if not dimensions in self.rose_space_domain_mappings:
            print("\t"+dimensions+" does not exist in space mappings?")
            import pdb; pdb.set_trace()

                
        #this_freq=rose_time_domain_mappings[freq]
        this_space=self.rose_space_domain_mappings[dimensions]

        #do these stash code(s) exist already in the existing STASH?

        #loop over all required stash variables for this cf variable
        plog("\n>> "+cf_var+" "+freq+" "+dims)
        #print(dims,this_space)
        for stash_code in req_stash:
           self.add_stash(stash_code,freq,dimensions)
           

           
    def rose_get_dom_level(self,dom_name):
        #returns the IOPL domain level index for dom_name in rose
        #get all the domain keys in rose
        rose_space_keys=[key for key in self.rose.keys() if 'umstash_domain' in key]
        for i in rose_space_keys:
            if self.rose[i]['dom_name']==dom_name:
                return(int(self.rose[i]['iopl']))
        print(dom_name+' not found?')
        import pdb; pdb.set_trace()


    def get_domain(self,stash_code,spatial_domain_cf):
        #checks to see if the domain defined by spatial_domain_cf matches what is required by stash_code
        #and returns the correct domain name, if possible
        #also now returns the spatial_domain_cf, as we sometimes have to modify this (mainly for check_output)
        spatial_domain=self.rose_space_domain_mappings[spatial_domain_cf]

        this_stash=self.stash_levels[stash_code]
        this_stash_level=self.rose_get_dom_level(spatial_domain)
        
        sc_level=this_stash['LevelT']
        sc_pseudo_level=this_stash['PseudT']

        
        #if sc_pseudo_level is NOT 0 then this stash_code uses pseudo levels - add to the spatial_domain_cf
        if sc_pseudo_level!=0:
      
            if not 'pseudo' in spatial_domain_cf:
                spatial_domain_cf+=' pseudo'
            
        
        #do we already have a mapping for this?
        if 'domains' in main_config:
            user_domains=main_config['domains']
            if stash_code in user_domains:

            
                #a mapping is defined in the config
                spatial_domain=user_domains[stash_code]
                plog("A domain mapping for "+stash_code+" was found in "+conf_file+": "+spatial_domain+" -will use this")
                return(spatial_domain,spatial_domain_cf)
            #there is a domains section defined in the config file
           

        #does this spatial_domain use the correct levels that this stash_code is output on?
        if this_stash_level!=sc_level:
            #No 
            if (this_stash_level==3) and (sc_level==2):
                #this request is on pressure levels - OK for theta and rho levels as they will be regridded
                plog("rho levels will be converted to pressure levels")
                import pdb; pdb.set_trace()
                #logging.info("rho levels will be converted to pressure levels")
            elif (this_stash_level==3) and (sc_level==1):
                #this request is on pressure levels - OK for theta and rho levels as they will be regridded
                plog("theta levels will be converted to pressure levels")
                import pdb; pdb.set_trace()

                #logging.info("theta levels will be converted to pressure levels")
            else:
                plog("Request is for "+stash_code+" using "+spatial_domain+" but this uses level="+str(this_stash_level)+" but "+str(stash_code)+" is output on level="+str(sc_level))
                #logging.info("Request is for "+stash_code+" using "+spatial_domain+" but this uses level="+str(this_stash_level)+" but "+str(stash_code)+" is output on level="+str(self.stash_levels[stash_code]))
                #import pdb; pdb.set_trace()
                if (this_stash_level==1) and (sc_level==2):
                    if spatial_domain==self.space_mappings['longitude latitude alevhalf time']:
                        plog("switching output from "+self.space_mappings['longitude latitude alevhalf time']+" to "+self.space_mappings['longitude latitude alevel time'])
                        #plog(ostr)
                        #logging.info(ostr)
                        spatial_domain_cf='longitude latitude alevel time'
                        spatial_domain=self.space_mappings[spatial_domain_cf]
                elif (this_stash_level==2) and (sc_level==1):
                    if spatial_domain==self.space_mappings['longitude latitude alevel time']:
                        plog("switching output from "+self.space_mappings['longitude latitude alevel time']+" to "+self.space_mappings['longitude latitude alevhalf time'])
                        #print(ostr)
                        #logging.info(ostr)
                        spatial_domain_cf='longitude latitude alevhalf time'
                        spatial_domain=self.space_mappings[spatial_domain_cf]
                elif (this_stash_level==1) and (sc_level==5):
                    #requested on RHO but this is a single level field - switch to 'DIAG'
                    plog("switching output from "+self.space_mappings['longitude latitude alevhalf time']+" to "+self.space_mappings['longitude latitude time'])
                    #print(ostr)
                    #logging.info(ostr)
                    spatial_domain_cf='longitude latitude time'
                    spatial_domain=self.space_mappings[spatial_domain_cf]
                elif (this_stash_level==2) and (sc_level==5):
                    #requested on RHO but this is a single level field - switch to 'DIAG'
                    plog("switching output from "+self.space_mappings['longitude latitude alevel time']+" to "+self.space_mappings['longitude latitude time'])  
                    #print(ostr)
                    #logging.info(ostr)
                    spatial_domain_cf='longitude latitude time'
                    spatial_domain=self.space_mappings[spatial_domain_cf]
                elif (this_stash_level==5) and (sc_level==6):
                    #this is probably a soil diagnostic that is being integrated over levels
                    #allow the level 6 (soil) levels to override
                    plog("Allow soil diagnostic to be output on all soil levels")
                    if "'DSOIL'" in self.rose_space_domain_mappings.values():
                        spatial_domain="'DSOIL'"
                        spatial_domain_cf+=' sdepth'
                        plog("Converted space domain to 'DSOIL'")

                    else:
                        print("Need DSOIL domain for this diagnostics, but not available???")
                        import pdb; pdb.set_trace()

                        
                else:
                    print("Not sure what to do here!")
                    import pdb; pdb.set_trace()


        if sc_pseudo_level!=0:
           #this field uses pseudo levels
           #which do we need?
           print(stash_code+" requires pseudo levels")

           
           rose_space_keys=[key for key in self.rose.keys() if 'umstash_domain' in key]
           pseudo_level_found=False
           for key in rose_space_keys:
              #compare the Pseudo level Type with the sc_pseudo_level
              if int(self.rose[key]['plt'])==sc_pseudo_level and int(self.rose[key]['iopl'])==sc_level:
                 if 'pslist' in self.rose[key]:
                    #OK this domain has a pslist
                    #create pslist of integers
                    this_dom=self.rose[key]['dom_name']
                    
                    pslist=[int(num) for num in self.rose[key]['pslist'].split(',')]
                    #Does this range of this pseudo level list for this domain match the range defined for the diagnostics in STASHMASTER?
                    if pslist[0]==this_stash['PseudF'] and pslist[-1]==this_stash['PseudL']:
                       #This domain should at least contain the required range defined for this diagnostic
                       #This may now always be correct - sometimes we might want a subset of pseudo levels - but how to specify this?
                       #eg only Plants or Trees in Surface tiles types
                       pseudo_level_found=True
                       break
                    else:
                       print("Found "+this_dom+" but the pseudolevel range "+str(pslist)+" does not match the defined range "+str(this_stash['PseudF'])+"-"+str(this_stash['PseudL'])+" for "+stash_code)
                       print("If you want to use "+this_dom+" for "+stash_code+" ("+self.stash_names[stash_code].strip()+") "+" add the line:\n"+stash_code+" = "+this_dom+"\nto "+main_config['user']['log_file'].strip("'")+" under a [domains] section")
                             
                    
           if pseudo_level_found:
              
              print(this_dom+" in ROSE uses the correct pseudo level and model levels "+key)
              plog("Using "+this_dom+" for "+stash_code)
              spatial_domain=this_dom
           else:
              print("Pseudo level range not found in ROSE domains")
              print("Checking CMIP reference")
            
              cmip6_space_keys=[key for key in self.cmip6.keys() if 'umstash_domain' in key]
              
              for key in cmip6_space_keys:
                 #compare the Pseudo level Type with the sc_pseudo_level
                 #also check that the IOPL is correct
                 if int(self.cmip6[key]['plt'])==sc_pseudo_level and  int(self.cmip6[key]['iopl'])==sc_level:
                    
                    if 'pslist' in self.cmip6[key]:
                       #OK this domain has a pslist
                       #create pslist of integers
                       this_dom=self.cmip6[key]['dom_name']
                       
                       pslist=[int(num) for num in self.cmip6[key]['pslist'].split(',')]
                       #Does this range of this pseudo level list for this domain match the range defined for the diagnostics in STASHMASTER?
                       if pslist[0]==this_stash['PseudF'] and pslist[-1]==this_stash['PseudL']:
                          #This domain should at least contain the required range defined for this diagnostic
                          #This may now always be correct - sometimes we might want a subset of pseudo levels - but how to specify this?
                          #eg only Plants or Trees in Surface tiles types
                          pseudo_level_found=True
                          break
                       else:
                          print("Found "+this_dom+" but the pseudolevel range "+str(pslist)+" does not match the defined range "+str(this_stash['PseudF'])+"-"+str(this_stash['PseudL'])+" for "+stash_code)
                    
                    break

              if pseudo_level_found:
                 print("Pseudo level "+str(sc_pseudo_level)+" found in CMIP6 in "+key)
                 import pdb; pdb.set_trace()
              else:
                 print("Pseudo level range  not found in CMIP6 Domains")
                 print("Need to define a domain usage in "+main_config['user']['log_file'].strip("'")+" under a [domains] section")
                 exit()
        return(spatial_domain,spatial_domain_cf)


    def get_time_domain(self,time_domain_cf,options,spatial_domain):
        '''
        find time domain for this cf time domain and the LBPROC in options
        '''

        if not 'lbproc' in options:
            lbproc='128'
            #is we can't find and LBPROC we assume this is just as time mean
        else:
            lbproc=options['lbproc']

        if not lbproc in self.lbproc_mappings:
            print("Hmmm.. "+lbproc+" not found!")
            import pdb; pdb.set_trace()
            #does rose or cmip6 contain a time domain that has this time operation
        this_time_domain_key=time_domain_cf+'_'+self.lbproc_mappings[lbproc]
        if this_time_domain_key in self.rose_time_domain_mappings:
           time_domain=self.rose_time_domain_mappings[this_time_domain_key]
        else:
           time_domain=''
        rose_time_keys=[key for key in self.rose.keys() if 'umstash_time' in key]
        cmip6_time_keys=[key for key in self.cmip6.keys() if 'umstash_time' in key]

        #pull in the part mapping for this time domain 
        time_filter=self.tim_dom[freq]
        #if ityp is not already defined in time_filterm then add the ityp for the method from the lbproc mappings
        if not 'ityp' in time_filter:
            time_filter['ityp']=self.lbproc_mappings[lbproc]
        rose_lbproc=[ x for x in rose_time_keys if is_subset(time_filter,um.rose[x])]
        time_usage_found=False
        #rose_lbproc=[ x for x in rose_time_keys if um.rose[x]['ityp']==self.lbproc_mappings[options['lbproc']]]
        if rose_lbproc:
           print("LBPROC found in ROSE")
           if len(rose_lbproc)>1:
              print("Hmm - we have more than one choice here!")
              import pdb; pdb.set_trace()
           else:
              time_domain=um.rose[rose_lbproc[0]]['tim_name']
              print("Switching to "+time_domain)
              return(time_domain)

        cmip6_lbproc=[ x for x in cmip6_time_keys if is_subset(time_filter,um.cmip6[x])]
        #cmip6_lbproc=[ x for x in cmip6_time_keys if um.cmip6[x]['ityp']==self.lbproc_mappings[options['lbproc']]]
        if cmip6_lbproc:
            print("LBPROC found in CMIP6")

            #COPY ACCROSS
            if len(cmip6_lbproc)>1:
                print("Hmm - we have more than one choice here!")
                if len(cmip6_lbproc)>2:
                    print("Too many!")
                    import pdb; pdb.set_trace()
                else:
                    if not is_equal_except(um.cmip6.items(cmip6_lbproc[0]),um.cmip6.items(cmip6_lbproc[1]),'tim_name'):
                        print("Two choices, but not equivalent!")
                        import pdb; pdb.set_trace()
                    print("2 choices, but these are equivalent time profiles")
                    print("We will use "+cmip6_lbproc[0])
             
                    
            time_domain=um.cmip6[cmip6_lbproc[0]]
            time_domain_name=time_domain['tim_name']
            #import pdb; pdb.set_trace()
            print("Copying "+time_domain_name+" to ROSE")
            self.copy_cmip6_tim_dom_to_rose(time_domain,spatial_domain)
            print("Switching to "+time_domain_name)

            return(time_domain_name)

        print("LBPROC not found anywhere!")
        import pdb; pdb.set_trace()

        return(time_domain) 

    def nc_check_stash(self,stash_code,time_domain_cf,spatial_domain_cf):
        print("Check output")
        #the split is because multiple occurrences of a stash code in a netcdf file will be appended with _2 _3 etc
        #we just want to check the stash code part
        #if 'pseudo' in spatial_domain_cf:
        #    import pdb; pdb.set_trace()
        
        matches=[key for key in nc_output if key.nc_get_variable().split('_')[0]==stash_code]
        if not matches:
            print(stash_code+" not found in NC output")
            spatial_domain_cf_list=sorted(spatial_domain_cf.split(' '))
        else:
            print(stash_code+" found in NC output")
            #if 'height' in spatial_domain_cf:
            #    import pdb; pdb.set_trace()

            #does this have the required domain?
            #We need to remap some of the requested dimension to what the model actually writes out
            replacements={'height10m':'',
                          'height2m':'',
                          'alevhalf':'model_level_number',
                          'alevel':'model_level_number',
                          'sdepth1':'depth',                          
                          'sdepth':'depth',
                          'typesi':''
                          }
            spatial_domain_cf_adjusted=spatial_domain_cf
            for torep in replacements:
                rep=replacements[torep]
                spatial_domain_cf_adjusted=spatial_domain_cf_adjusted.replace(torep,rep)

            
            #strip out any empty strings following replacement                                                             

            #.replace('height10m','height').replace('height2m','height').replace('alevhalf','model_level_number').replace('alevel','model_level_number')
            spatial_domain_cf_list=sorted(spatial_domain_cf_adjusted.split(' '))
            #spatial_domain_cf_list=sorted(spatial_domain_cf.split(' '))

            #strip out any empty strings following replacement                                                             
            spatial_domain_cf_list=[item for item in spatial_domain_cf_list if item!='']


            for this_match in matches:
                nc_domain=sorted([item.identity() for item in this_match.coords().values()])
                #nc_domain=sorted([item.standard_name for item in this_match.coords().values()])
                if 'air_pressure' in nc_domain:
                    # if the domain contains air pressure - let's guess what the original plev was!
                    new_name='plev'+str(this_match.coord('air_pressure').size)
                    nc_domain=sorted([item if item!='air_pressure' else new_name for item in nc_domain])
                if 'long_name=Land and Vegetation Surface types' in nc_domain:
                    # if the domain contains veg and surface types -this is a pseudo level
                    new_name='pseudo'
                    nc_domain=sorted([item if item!='long_name=Land and Vegetation Surface types' else new_name for item in nc_domain])

                if 'height' in nc_domain:
                    #domain contains a height coordinate
                    if this_match.coord('height').size==1:
                        #this is a single level - hence we can ignore here
                        nc_domain=sorted([item for item in nc_domain if item!='height'])
                #if 'model_level_number' in nc_domain:
                #    nc_domain=sorted([item if item!='model_level_number' else 'alevel' for item in nc_domain])
                #    nc_domain=sorted([item if item!='model_level_number' else 'alevhalf' for item in nc_domain])
                if spatial_domain_cf_list==nc_domain:
                    #spatial domains matc
                    this_time_domain=this_match.properties()['name'].split('_')[1].replace('6hpt','6hrPt')
                    if time_domain_cf == this_time_domain:
                        #time domains match
                        print("Time and spatial domains match")
                        um.nc_found.append([stash_code,time_domain_cf,spatial_domain_cf])
                        return(True)
        #if we get to here, there were no matches!
        print("Couldn't find "+stash_code+" in NC for "+' '.join(spatial_domain_cf_list)+" at "+time_domain_cf)
        self.nc_missing.append([stash_code,time_domain_cf,spatial_domain_cf])
        #import pdb; pdb.set_trace()

        return(False)
     
    def add_stash(self,stash_code0,time_domain_cf,spatial_domain_cf):
        #checks to see if the stash_code (e.g m01s03i236) exists in the rose config object
        #and that this stash code is output using the correct time and spatial domains
        #if not, adds the stash code and any required domains



        pattern = r'm(\d{2})s(\d{2})i(\d{3})\[(.*?)\]'
        if re.match(pattern,stash_code0):
            model,isec,item,options=re.match(pattern,stash_code0).groups()
            options=dict(pair.split('=') for pair in options.split(','))
        else:
            pattern = r'm(\d{2})s(\d{2})i(\d{3})'
            model,isec,item=re.match(pattern,stash_code0).groups()
            options={}
        stash_code='m'+model+'s'+isec+'i'+item


        
        
        #rearrnge options as a dict

        #default time domains
        

        spatial_domain,spatial_domain_cf=self.get_domain(stash_code,spatial_domain_cf)

        if check_output:
            self.nc_check_stash(stash_code,time_domain_cf,spatial_domain_cf)
            #if this returns true - the the stash code exists with this time and space domain in the netcdf
            return()

        plog("Add "+stash_code)
        
        time_domain=self.get_time_domain(time_domain_cf,options,spatial_domain)
        #import pdb; pdb.set_trace()

        

        if 'blev' in options:
           print("BLEV")
           blev_lc= options['blev'].lower()
           if blev_lc=='0.05':
              #this maps onto the first soil level
              blev_lc='sdepth1'
           if blev_lc in dims:
                print(stash_code+" required "+  blev_lc +" - and this DOES exist in "+dims)
           else:
                print(stash_code+" required "+  blev_lc +" - but this DOES NOT exist in "+dims)

                groups=re.match(r'.*?\s(plev.*?)\s.*?',dims)
                if not groups:
                   print("Not sure what to do now!")
                   import pdb; pdb.set_trace()

                else:
                   print("However, the levels request is for "+groups[1]+" so we will continue with this")
                
        
        
        
        this_stash=self.stash_levels[stash_code]
            
        #pattern = r'm(\d{2})s(\d{2})i(\d{3})'
        #model,isec,item=re.match(pattern,stash_code).groups()

        
        this_stash_level=self.rose_get_dom_level(spatial_domain)
        sc_level=this_stash['LevelT']
        sc_pseudo_level=this_stash['PseudT']

               
        if model!="01":
            print("Unknown model in stashcode?")
            print(stash_code)
            import pdb; pdb.set_trace()


        stash_found=False
        time_found=False
        space_found=False
        #extract all the umstash requests from rose 
        stash_c= [key for key in self.rose.sections() if re.search('umstash_streq', key)]
        #loop over all these and see if one matches the stash id we require at the same freq output and domain
        for req in stash_c:
            this_stash=self.rose[req]
            #if the isec and item match and this is the atmospher (model=01)
            if (this_stash['isec'].zfill(2)==isec) and (this_stash['item'].zfill(3)==item) and model=='01':

                #if 'm01s03i460' in stash_code:
                #   import pdb; pdb.set_trace()
                this_space=this_stash['dom_name']
                this_time=this_stash['tim_name']
                #does this stash have the required space and time domain?
                
                if (this_time==time_domain) and (this_space==spatial_domain):
                    stash_found=True
                    plog(stash_code+" already exists in "+self.rose_stash+" with "+time_domain+" and "+spatial_domain+" so no need to add anything")
                    #logging.info(stash_code+" already exists at "+time_domain+" and "+spatial_domain)
                    
                    #if check_output:
                    #    import pdb; pdb.set_trace()
                    #   if self.nc_check_stash(stash_code,time_domain_cf,spatial_domain_cf):
                    #        #if this returns true - the the stash code exists with this time and space domain in the netcdf
                    #        return()
                            
                            
                    break

        #if the stash was found, it already exists, so we don't need to do anything - otherwise..
        if not stash_found:
            #stash code not found in rose with correct time and space domain:
            print("Need to add stash code "+stash_code)
            

            #for USAGE - we will pick the FIRST usage in the list from use_matrix[time][space]
            #Does the time_domain exist in the use_matrix?
            #if 'm01s03i236' in stash_code:
            #   import pdb; pdb.set_trace()

            if not time_domain in self.use_matrix:
                print(time_domain+" doesn't exist in the use matrix?")
                #can we find this time domain in CMIP6?
                import pdb; pdb.set_trace()

                #what about the spatial_domain?
            elif not spatial_domain in self.use_matrix[time_domain]:
                plog(spatial_domain+" doesn't exist in the use matrix?")
                #OK so this domain is not associated with an existing USAGE
                #is there a DEFAULT for this time_domain?


                if time_domain_cf in self.default_usage:
                    #there is!
                    usage=self.default_usage[time_domain_cf]
                    plog("using the default usage for "+time_domain+" of "+usage)
                    #check that this default usage actually exists in rose_conf - if not - take action!
                    self.check_use_already_exists(usage)
                    
                    #Can't assume that all STASH have "UPD" defined as ann output stream
                    #

                    #print(oft)
                    #logging.info(oft)
                else:
                    #Ok, so no default usage and no spatial_domain entry in the use_matrux
                    #so we will pick the first one that is associated with this time_domain
                    u1=[]
                    for i in self.use_matrix[time_domain].values():
                        u1.extend(i)
                    u2=list(set(u1))
                    u2.sort()
                    usage=u2[0]
                    plog("No default usage for "+time_domain+" picking first: "+usage)

            else:
                usage=self.use_matrix[time_domain][spatial_domain][0]

            plog("Will use "+usage+" as the usage")

            #print(oft)
            #logging.info(oft)
            isec1=isec.lstrip('0') if isec!='00' else '0'
            item1=item.lstrip('0') if item!='000' else '0'




            # Convert the UUID to an 8-digit hexadecimal representation
            #hex_uuid = format(int(uuid.uuid4().hex[:8], 16), 'x')
            
            #Added ens_name here, but maybe this just needs to be added if we are writing to XIOS rather than UM STASH?
            new_stash={'dom_name':spatial_domain,
                       'ens_name':"''",
                       'isec':isec1,
                       'item':item1,
                       'package':"'EXTRA'",
                       'tim_name':time_domain,
                       'use_name':usage
                       }

            #computes the correct hash uuid for this stash
            hex_uuid=self.get_uuid_hash(new_stash)
            #namelist_name="[!namelist:umstash_streq("+isec+item+"_"+hex_uuid+")]"
            namelist_name="namelist:umstash_streq("+isec+item+"_"+hex_uuid+")"
            
            self.rose[namelist_name]=new_stash
            plog("Added new stash entry for "+stash_code+" to ROSE using "+spatial_domain+" and "+time_domain)
            self.added.append(stash_code)
            #logging.info('Added '+stash_code+' using '+time_domain+' '+spatial_domain+' '+usage)
        return(stash_found)

        
    def add_stash_old(self,stash_code,time_domain_cf,spatial_domain_cf):
        #checks to see if the stash_code (e.g m01s03i236) exists in the rose config object
        #and that this stash code is output using the correct time and spatial domains
        #if not, adds the stash code and any required domains

        spatial_domain=self.get_domain(stash_code,spatial_domain_cf)
        
       
        time_domain=self.rose_time_domain_mappings[time_domain_cf]
        spatial_domain=self.rose_space_domain_mappings[spatial_domain_cf]

        this_stash=self.stash_levels[stash_code]
        #extract all stash requests from rose
        pattern = r'm(\d{2})s(\d{2})i(\d{3})'
        model,isec,item=re.match(pattern,stash_code).groups()
        this_stash_level=self.rose_get_dom_level(spatial_domain)
        sc_level=this_stash['LevelT']
        sc_pseudo_level=this_stash['PseudT']

        
        
        #does this spatial_domain use the correct levels that this stash_code is output on?
        if this_stash_level!=sc_level:
            if (this_stash_level==3) and (sc_level==2):
                #this request is on pressure levels - OK for theta and rho levels as they will be regridded
                plog("rho levels will be converted to pressure levels")
                #logging.info("rho levels will be converted to pressure levels")
            elif (this_stash_level==3) and (sc_level==1):
                #this request is on pressure levels - OK for theta and rho levels as they will be regridded
                plog("theta levels will be converted to pressure levels")
                #logging.info("theta levels will be converted to pressure levels")
            else:
                plog("Request is for "+stash_code+" using "+spatial_domain+" but this uses level="+str(this_stash_level)+" but "+str(stash_code)+" is output on level="+str(sc_level))
                #logging.info("Request is for "+stash_code+" using "+spatial_domain+" but this uses level="+str(this_stash_level)+" but "+str(stash_code)+" is output on level="+str(self.stash_levels[stash_code]))
                if (this_stash_level==1) and (sc_level==2):
                    if spatial_domain==self.space_mappings['longitude latitude alevhalf time']:
                        plog("switching output from "+self.space_mappings['longitude latitude alevhalf time']+" to "+self.space_mappings['longitude latitude alevel time'])
                        #plog(ostr)
                        #logging.info(ostr)
                        spatial_domain=self.space_mappings['longitude latitude alevel time']
                elif (this_stash_level==2) and (sc_level==1):
                    if spatial_domain==self.space_mappings['longitude latitude alevel time']:
                        plog("switching output from "+self.space_mappings['longitude latitude alevel time']+" to "+self.space_mappings['longitude latitude alevhalf time'])
                        #print(ostr)
                        #logging.info(ostr)
                        spatial_domain=self.space_mappings['longitude latitude alevhalf time']
                elif (this_stash_level==1) and (sc_level==5):
                    #requested on RHO but this is a single level field - switch to 'DIAG'
                    plog("switching output from "+self.space_mappings['longitude latitude alevhalf time']+" to "+self.space_mappings['longitude latitude time'])
                    #print(ostr)
                    #logging.info(ostr)

                    spatial_domain=self.space_mappings['longitude latitude time']
                elif (this_stash_level==2) and (sc_level==5):
                    #requested on RHO but this is a single level field - switch to 'DIAG'
                    plog("switching output from "+self.space_mappings['longitude latitude alevel time']+" to "+self.space_mappings['longitude latitude time'])
                    #print(ostr)
                    #logging.info(ostr)

                    spatial_domain=self.space_mappings['longitude latitude time']
                elif (this_stash_level==5) and (sc_level==6):
                    #this is probably a soil diagnostic that is being integrated over levels
                    #allow the level 6 (soil) levels to override
                    plog("Allow soil diagnostic to be output on all soil levels")
                    if "'DSOIL'" in self.rose_space_domain_mappings.values():
                        spatial_domain="'DSOIL'"
                        plog("Converted space domain to 'DSOIL'")
                    else:
                        print("Need DSOIL domain for this diagnostics, but not available???")
                        import pdb; pdb.set_trace()

                        
                else:
                    print("Not sure what to do here!")
                    import pdb; pdb.set_trace()


        if sc_pseudo_level!=0:
           #this field uses pseudo levels
           #which do we need?
           print(stash_code+" requires pseudo levels")

           
           rose_space_keys=[key for key in self.rose.keys() if 'umstash_domain' in key]
           pseudo_level_found=False
           for key in rose_space_keys:
              #compare the Pseudo level Type with the sc_pseudo_level
              if int(self.rose[key]['plt'])==sc_pseudo_level and int(self.rose[key]['iopl'])==sc_level:
                 if 'pslist' in self.rose[key]:
                    #OK this domain has a pslist
                    #create pslist of integers
                    this_dom=self.rose[key]['dom_name']
                    
                    pslist=[int(num) for num in self.rose[key]['pslist'].split(',')]
                    #Does this range of this pseudo level list for this domain match the range defined for the diagnostics in STASHMASTER?
                    if pslist[0]==this_stash['PseudF'] and pslist[-1]==this_stash['PseudL']:
                       #This domain should at least contain the required range defined for this diagnostic
                       #This may now always be correct - sometimes we might want a subset of pseudo levels - but how to specify this?
                       #eg only Plants or Trees in Surface tiles types
                       pseudo_level_found=True
                       break
                    else:
                       print("Found "+this_dom+" but the pseudolevel range "+str(pslist)+" does not match the defined range "+str(this_stash['PseudF'])+"-"+str(this_stash['PseudL'])+" for "+stash_code)
                       print("If you want to use "+this_dom+" for "+stash_code+" add the line:\n"+stash_code+" = "+this_dom+"\nto "+main_config['user']['log_file'].strip("'")+" under a [domains] section")
                             
                    
           if pseudo_level_found:
              
              print(this_dom+" in ROSE uses the correct pseudo level and model levels "+key)
              plog("Using "+this_dom+" for "+stash_code)
              spatial_domain=this_dom
           else:
              print("Pseudo level range not found in ROSE domains")
              print("Checking CMIP reference")
            
              cmip6_space_keys=[key for key in self.cmip6.keys() if 'umstash_domain' in key]
              
              for key in cmip6_space_keys:
                 #compare the Pseudo level Type with the sc_pseudo_level
                 #also check that the IOPL is correct
                 if int(self.cmip6[key]['plt'])==sc_pseudo_level and  int(self.cmip6[key]['iopl'])==sc_level:
                    
                    if 'pslist' in self.cmip6[key]:
                       #OK this domain has a pslist
                       #create pslist of integers
                       this_dom=self.cmip6[key]['dom_name']
                       
                       pslist=[int(num) for num in self.cmip6[key]['pslist'].split(',')]
                       #Does this range of this pseudo level list for this domain match the range defined for the diagnostics in STASHMASTER?
                       if pslist[0]==this_stash['PseudF'] and pslist[-1]==this_stash['PseudL']:
                          #This domain should at least contain the required range defined for this diagnostic
                          #This may now always be correct - sometimes we might want a subset of pseudo levels - but how to specify this?
                          #eg only Plants or Trees in Surface tiles types
                          pseudo_level_found=True
                          break
                       else:
                          print("Found "+this_dom+" but the pseudolevel range "+str(pslist)+" does not match the defined range "+str(this_stash['PseudF'])+"-"+str(this_stash['PseudL'])+" for "+stash_code)
                    
                    break

              if pseudo_level_found:
                 print("Pseudo level "+str(sc_pseudo_level)+" found in CMIP6 in "+key)
                 import pdb; pdb.set_trace()
              else:
                 print("Pseudo level range  not found in CMIP6 Domains")
                 print("Need to define a domain usage in "+main_config['user']['log_file'].strip("'")+" under a [domains] section")
                 exit()
        
        if model!="01":
            print("Unknown model in stashcode?")
            print(stash_code)
            import pdb; pdb.set_trace()


        stash_found=False
        time_found=False
        space_found=False
        #extract all the umstash requests from rose 
        stash_c= [key for key in self.rose.sections() if re.search('umstash_streq', key)]
        #loop over all these and see if one matches the stash id we require at the same freq output and domain
        for req in stash_c:
            this_stash=self.rose[req]
            #if the isec and item match and this is the atmospher (model=01)
            if (this_stash['isec'].zfill(2)==isec) and (this_stash['item'].zfill(3)==item) and model=='01':

                this_space=this_stash['dom_name']
                this_time=this_stash['tim_name']
                #does this stash have the required space and time domain?
     
                if (this_time==time_domain) and (this_space==spatial_domain):
                    stash_found=True
                    plog(stash_code+" found in ROSE with "+time_domain+" and "+spatial_domain)
                    #logging.info(stash_code+" already exists at "+time_domain+" and "+spatial_domain)
                    break

        #if the stash was found, it already exists, so we don't need to do anything - otherwise..
        if not stash_found:
            #stash code not found in rose with correct time and space domain:
            print("Need to add stash code "+stash_code)
           

            #for USAGE - we will pick the FIRST usage in the list from use_matrix[time][space]
            #Does the time_domain exist in the use_matrix?
            if not time_domain in self.use_matrix:
                print(time_domain+" doesn't exist in the use matrix?")
                #can we find this time domain in CMIP6?
                import pdb; pdb.set_trace()

                #what about the spatial_domain?
            elif not spatial_domain in self.use_matrix[time_domain]:
                plog(spatial_domain+" doesn't exist in the use matrix?")
                #OK so this domain is not associated with an existing USAGE
                #is there a DEFAULT for this time_domain?


                if time_domain_cf in self.default_usage:
                    #there is!
                    usage=self.default_usage[time_domain_cf]
                    plog("using the default usage for "+time_domain+" of "+usage)
                    #check that this default usage actually exists in rose_conf - if not - take action!
                    self.check_use_already_exists(usage)
                    
                    #Can't assume that all STASH have "UPD" defined as ann output stream
                    #

                    #print(oft)
                    #logging.info(oft)
                else:
                    #Ok, so no default usage and no spatial_domain entry in the use_matrux
                    #so we will pick the first one that is associated with this time_domain
                    u1=[]
                    for i in self.use_matrix[time_domain].values():
                        u1.extend(i)
                    u2=list(set(u1))
                    u2.sort()
                    usage=u2[0]
                    plog("No default usage for "+time_domain+" picking first: "+usage)

            else:
                usage=self.use_matrix[time_domain][spatial_domain][0]

            plog("Will use "+usage+" as the usage")

            #print(oft)
            #logging.info(oft)
            isec1=isec.lstrip('0') if isec!='00' else '0'
            item1=item.lstrip('0') if item!='000' else '0'




            # Convert the UUID to an 8-digit hexadecimal representation
            #hex_uuid = format(int(uuid.uuid4().hex[:8], 16), 'x')
            
            #Added ens_name here, but maybe this just needs to be added if we are writing to XIOS rather than UM STASH?
            new_stash={'dom_name':spatial_domain,
                       'ens_name':"''",
                       'isec':isec1,
                       'item':item1,
                       'package':"'EXTRA'",
                       'tim_name':time_domain,
                       'use_name':usage
                       }

            #computes the correct hash uuid for this stash
            hex_uuid=self.get_uuid_hash(new_stash)
            #namelist_name="[!namelist:umstash_streq("+isec+item+"_"+hex_uuid+")]"
            namelist_name="namelist:umstash_streq("+isec+item+"_"+hex_uuid+")"
            
            self.rose[namelist_name]=new_stash
            plog("Added new stash entry for "+stash_code+" to ROSE")
            self.added.append(stash_code)
            #logging.info('Added '+stash_code+' using '+time_domain+' '+spatial_domain+' '+usage)
        return(stash_found)

    def write(self,rose_outfile):
        if self.rose_header!='':
            rose_out=open(rose_outfile,'w')
            rose_out.write(self.rose_header+"\n\n")
            rose_out.close()

            with open(rose_outfile, 'a') as rose_out:
                self.rose.write(rose_out)
        else:
            with open(rose_outfile, 'w') as rose_out:
                self.rose.write(rose_out)
            
        print("Written "+rose_outfile)



class CICE:


    def __init__(self):


        self.freq_map={'1d':'day','1m':'mon'}
        #self.read_rose_app_conf(file+'/'+ice_conf)
        self.rose_cice=main_config['user']['job_path']+'app/nemo_cice/rose-app.conf'
        self.cice_diagnostics_file=main_config['main']['cice_diags']
        self.rose,self.rose_header=self.read_rose_app_conf(self.rose_cice)

        self.nc_found=[] #list of UM diagnostics found in the NC file during --check_output
        self.nc_missing=[] #list of UM diagnostics missing in the NC file during --check_output

        self.missing=[] # list of diagnostics we failed to add!
        self.added=[]#list of added diagnostics

        self.cice_diagnostics=self.read_cice_diagnostics(self.cice_diagnostics_file)
        

    def read_cice_diagnostics(self,file):
       #file='ice_history_shared.F90'

       namelist=False
       with open(file, 'r') as infile:
          lines=infile.readlines()
          cice_diagnostics=[]
          for line in lines:
             if namelist:
                cice1=line.replace(' ','').replace('\t','').replace('\n','').split(',')
                cice2=[i for i in cice1 if not '&' in i and not '!' in i]
                cice_diagnostics.extend(cice2)
                if not '&' in line:
                   namelist=False
             if 'icefields_nml' in line and 'namelist' in line:
                print(line)
                namelist=True
       return(cice_diagnostics)

  
        
    def read_rose_app_conf(self,file):
        config = configparser.ConfigParser()
        #turn of the lowercaseization of keys
        config.optionxform = lambda option: option
        configFilePath = file
        #Need to do this as the header in rose-app.conf files doesn't have a [DEFAULT] section heading
        header_flag=True
        with open(configFilePath) as stream:
            header,config_string=stream.read().split('\n',1)

            config.read_string(config_string)
            #config.read_string("[DUMMY]\n" + stream.read())
        return(config,header)


    def contains_operators(self,input_string):
       operators = set("*+-/")

       # Check if any of the special characters are in the string
       return any(i in operators for i in input_string)


    def nc_check_ice(self,fdiag,freq,dims):
        print("Check Ice output")
        #CICE doesn't have sea ice types, or classes, so remove this dimension
        if 'present' in fdiag:
            fdiag='f_ice_present'

        dims=dims.replace('typesi','')
        diag=fdiag.replace('f_','')
        matches=[key for key in nc_output if key.nc_get_variable()==diag]
        if matches:
            print(diag+" found in NC output")
            #import pdb; pdb.set_trace()
            #does this have the required domain?
            spatial_domain_cf_list=sorted(dims.split())
            ##HERE
            ## SPATIAL domain doesn't map perfectly as ICE is IJ not long lat!
            for this_match in matches:
                nc_domain=[item.nc_get_variable() for item in this_match.coords().values()]
                nc_domain=sorted([item.replace('TLON','longitude').replace('TLAT','latitude') for item in nc_domain])
                nc_domain=sorted([item.replace('ULON','longitude').replace('ULAT','latitude') for item in nc_domain])
                nc_domain=sorted([item.replace('VLON','longitude').replace('VLAT','latitude') for item in nc_domain])
                #nc_domain=sorted([item.standard_name for item in this_match.coords().values()])
                if spatial_domain_cf_list==nc_domain:
                    #spatial domains match!
                    ##THIS DOESN'T work for CICE
                    #try and get the frequency of this variable from the original file name - should be 1d or 1m
                    match_freq=this_match.get_filenames().pop().split('_')[-2]
                    if freq == self.freq_map[match_freq]:
                        #time domains match
                        print("Time and spatial domains match")
                        self.nc_found.append([diag,freq,dims])
                        return(True)

        else:
            print(diag+" not found in NC output")

        print("Couldn't find "+diag+" in NC for "+dims+" at "+freq)
        self.nc_missing.append([diag,freq,dims])
        return(False)
    
    def addIceDiag(self,fdiag,freq,dims):
    
        print("")

        print("-------------------------")
        #plog("Ice: "+diag)
        

        #check frequency is valid
        freq_map={'mon':"'m'",'day':"'d'",'hour':"'h'",'timestep':"'1'",'year':"'y'"}
        if not freq in freq_map:
            print("unknown output frequency"+freq+' '+fdiag)
            import pdb; pdb.set_trace()
        this_freq=freq_map[freq]

        if check_output:
            if self.nc_check_ice(fdiag,freq,dims):
                return()
                import pdb; pdb.set_trace()

        
        #if 'ice_present' in fdiag:
        #    import pdb; pdb.set_trace()
        #
        #           #mappings in common_mappings.cfg is wrong for sitimefrac !
        #ice_present in the confing by ice_history_shared.F90 is f_icepresent
        #           fdiag='f_icepresent'
        #        else:
        #           fdiag='f_'+this_expression

                    
        this_section=self.rose['namelist:icefields_nml']
        if fdiag in this_section:
           #this IS a valid CICE diagnostic

           diag_freq=this_section[fdiag]

           #this_freq = something line "'m'" so need to remove ' for comparison to work
           if not this_freq.strip("'") in diag_freq:
              plog(fdiag+" not currently output at "+freq+" - adding..")
              #if diag_freq is just 'x' we replace with this_freq
              if "x" in diag_freq:
                 diag_freq=this_freq
                 self.added.append(fdiag)
              else:
                 #otherwise, we need to add the new freq to the string (making sure it is surrounded by '  ')
                 diag_freq=(diag_freq+this_freq).replace("\'\'","")
                 this_section[fdiag]=diag_freq
                 self.added.append(fdiag)
                 #overwrite this_section in the ice dict
                 #we don't need to do this, as this_section already is a reference to self.rose['namelist:icefields_nml']
                 #self.rose['namelist:icefields_nml']=this_section
              #is diag_freq already set in the histfreq section of setup_nml?
              setup=self.rose['namelist:setup_nml']
              histfreq=setup['histfreq'].split(',')
              histfreq_n=setup['histfreq_n'].split(',')
              #is the output frequency for this diagnostic already present in the histfreq variable?

              if not this_freq in histfreq:
                 #no - so we need to add it to the next available slot
                 plog(this_freq+' not in current histfreq: '+','.join(histfreq))
                 found=False
                 for i in range(len(histfreq)):
                    if histfreq[i]=="'x'":
                       histfreq[i]=this_freq
                       #we also need to add how often this frequency needs to be output
                       histfreq_n[i]='1'
                       found=True
                       #if we are outputting daily, then the first element in histfreq should be "'d'"
                       #this is so that postproc picks up the daily files and concantenates them into a single monthly file of daily means
                       if "'d'" in histfreq:
                           d_pos=histfreq.index("'d'")
                           histfreq=["'d'"]+histfreq[0:d_pos]+histfreq[d_pos+1:]
                           histfreq_n=[histfreq_n[d_pos]]+histfreq_n[0:d_pos]+histfreq_n[d_pos+1:]

                       break
                 if not found:
                    print('no space for addditional CICE diag frequency!')
                    import pdb; pdb.set_trace()
                 setup['histfreq']=','.join(histfreq)
                 setup['histfreq_n']=','.join(histfreq_n)
                 #also don;t need to do this!
                 #self.rose['namelist:setup_nml']=setup

           else:
              plog(fdiag+" is already output at "+freq)
              

        else:
          
           if fdiag in self.cice_diagnostics:
              plog(fdiag+" not in the current CICE diag namelist - adding")
              import pdb; pdb.set_trace()

           else:
              plog(bold(fdiag+" not found - SKIPPING"))
              self.missing.append(diag)
       

        return()

    def addIceDiag_old2(self,line):
        diag=line['variable']
        freq=line['time']
        dims=line['space']

        print("")

        print("-------------------------")
        #plog("Ice: "+diag)

        #check frequency is valid
        freq_map={'mon':"'m'",'day':"'d'",'hour':"'h'",'timestep':"'1'",'year':"'y'"}
        if not freq in freq_map:
            print("unknown output frequency"+freq+' '+diag)
            import pdb; pdb.set_trace()
        this_freq=freq_map[freq]

        
        #is this a valid CICE diagnostics?
        


        #is this diag in the cf mapping tables?
        if not diag in um.cf_to_stash:
           #No - so let's guess it is like the CICE diags
           fdiag='f_'+diag
        else:
           #OK attempt to extract
           this_map=um.cf_to_stash[diag]
           this_expression=this_map['expression']
           pattern = r'm\d{2}s\d{2}i\d{3}\[.*?\]'
           matches = re.findall(pattern, this_expression)
           if len(matches)>0:
              #Aha - this IS a sea ice diagnostics, but the mapping files use atmosphere diagnostics to build in
              #So call the UM diagnostics function, and then return
              print("We need to add UM diagnostics for this "+diag+" CICE diagnostic!")
              um.add_cf_diagnostic(line)
              print("Done")
              return()
           else:
              if self.contains_operators(this_expression):
                 #This is some expression involving CONSTANTS (upper case) diagnostics (lower case), operators (*+.-) and brackets
                 #e.g. ((meltt + meltb + meltl - congel - frazil) * ICE_DENSITY + melts * SNOW_DENSITY) / (100. * SECONDS_IN_DAY * aice) * REF_SALINITY / 1000.
                 # we need to extract JUST the diagnostics
                 translation_table = str.maketrans("", "", "*+-/().0123456789_ABCDEFGHIJKLMNOPQRSTUVWXYZ")
                 #strip out all operators, numbers and upper case characters
                 sub_diag_string=this_expression.translate(translation_table)
                 #split into a list of diagnostics
                 sub_diags=list(filter(None, sub_diag_string.split()))
                 if len(sub_diags)==0:
                    print("something went wrong - mapping expression failed!")
                    import pdb; pdb.set_trace()
                 elif len(sub_diags)>1:
                    print(diag+" maps to more than one ice diagnostic - loop over these")
                    for sub_diag in sub_diags:
                       #call this method again
                       #take a copy of line
                       new_line=dict(line)
                       #swap variable
                       new_line['variable']=sub_diag
                       #call addIceDiag again
                       print("Adding "+sub_diag)
                       self.addIceDiag(new_line)                       
                    print("Done")
                    return()
                 else:
                    fdiag='f_'+sub_diags[0]
              else:
                 if this_expression=="ice_present":
                    #mappings in common_mappings.cfg is wrong for sitimefrac !
                    #ice_present in the confing by ice_history_shared.F90 is f_icepresent
                    fdiag='f_icepresent'
                 else:
                    fdiag='f_'+this_expression
           

        this_section=self.rose['namelist:icefields_nml']
        if fdiag in this_section:
           #this IS a valid CICE diagnostic

           diag_freq=this_section[fdiag]

           #this_freq = something line "'m'" so need to remove ' for comparison to work
           if not this_freq.strip("'") in diag_freq:
              plog(diag+" not currently output at "+freq+" - adding..")
              #if diag_freq is just 'x' we replace with this_freq
              if "x" in diag_freq:
                 diag_freq=this_freq
                 self.added.append(diag)
              else:
                 #otherwise, we need to add the new freq to the string (making sure it is surrounded by '  ')
                 diag_freq=(diag_freq+this_freq).replace("\'\'","")
                 this_section[fdiag]=diag_freq
                 self.added.append(diag)
                 #overwrite this_section in the ice dict
                 #we don't need to do this, as this_section already is a reference to self.rose['namelist:icefields_nml']
                 #self.rose['namelist:icefields_nml']=this_section
              #is diag_freq already set in the histfreq section of setup_nml?
              setup=self.rose['namelist:setup_nml']
              histfreq=setup['histfreq'].split(',')
              histfreq_n=setup['histfreq_n'].split(',')
              #is the output frequency for this diagnostic already present in the histfreq variable?

              if not this_freq in histfreq:
                 #no - so we need to add it to the next available slot
                 plog(this_freq+' not in current histfreq: '+','.join(histfreq))
                 found=False
                 for i in range(len(histfreq)):
                    if histfreq[i]=="'x'":
                       histfreq[i]=this_freq
                       #we also need to add how often this frequency needs to be output
                       histfreq_n[i]='1'
                       found=True
                       break
                 if not found:
                    print('no space for addditional CICE diag frequency!')
                    import pdb; pdb.set_trace()
                 setup['histfreq']=','.join(histfreq)
                 setup['histfreq_n']=','.join(histfreq_n)
                 #also don;t need to do this!
                 #self.rose['namelist:setup_nml']=setup

           else:
              plog(diag+" is already output at "+freq)
              

        else:
          
           if fdiag in self.cice_diagnostics:
              plog(diag+" not in the current CICE diag namelist - adding")
              import pdb; pdb.set_trace()

           else:
              plog(bold(diag+" not found - SKIPPING"))
              self.missing.append(diag)
       

        return()
     
    def addIceDiag_old(self,line):
        diag=line['variable']
        freq=line['time']
        dims=line['space']

        print("")

        print("-------------------------")
        #plog("Ice: "+diag)
        

        if ' ' in freq:
            import pdb; pdb.set_trace()

        freq_map={'mon':"'m'",'day':"'d'",'hour':"'h'",'timestep':"'1'",'year':"'y'"}
        if not freq in freq_map:
            print("unknown output frequency"+freq+' '+diag)
            import pdb; pdb.set_trace()
        this_freq=freq_map[freq]
        setup=self.rose['namelist:setup_nml']
        histfreq=setup['histfreq'].split(',')
        histfreq_n=setup['histfreq_n'].split(',')
        #is the output frequency for this diagnostic already present in the histfreq variable?
        if not this_freq in histfreq:
            #no - so we need to add it to the next available slot
            plog(this_freq+' not in current histfreq: '+','.join(histfreq))
            found=False
            for i in range(len(histfreq)):
                if histfreq[i]=="'x'":
                    histfreq[i]=this_freq
                    #we also need to add how often this frequency needs to be output
                    histfreq_n[i]='1'
                    found=True
                    break
            if not found:
                print('no space for addditional CICE diag frequency!')
                import pdb; pdb.set_trace()
            #write modified histfreq and histfreq_n back to ice
            setup['histfreq']=','.join(histfreq)
            setup['histfreq_n']=','.join(histfreq_n)
            self.rose['namelist:setup_nml']=setup

        import pdb; pdb.set_trace()
        this_map=um.cf_to_stash[diag]
        fdiag='f_'+diag
        this_section=self.rose['namelist:icefields_nml']
        if fdiag in this_section:
            diag_freq=this_section[fdiag]
            if not this_freq in diag_freq:
                plog(diag+" not currently output at "+freq+" - adding..")
                #if diag_freq is just 'x' we replace with this_freq
                if "x" in diag_freq:
                    diag_freq=this_freq
                else:
                    #otherwise, we need to add the new freq to the string (making sure it is surrounded by '  ')
                    diag_freq=(diag_freq+this_freq).replace("\'\'","")
                this_section[fdiag]=diag_freq
                #overwrite this_section in the ice dict
                self.rose['namelist:icefields_nml']=this_section

            else:
                plog(diag+" is already output at "+freq)


        else:
            plog(bold(diag+" not found - SKIPPING"))
            self.missing.append(diag)



        return()

    
    def write(self,rose_outfile):
        rose_out=open(rose_outfile,'w')
        rose_out.write(self.rose_header+"\n\n")
        rose_out.close()

        with open(rose_outfile, 'a') as rose_out:
            self.rose.write(rose_out)

        print("Written "+rose_outfile)


        
    
    pass


################### OCEAN

class Nemo:
    #NEMO diagnostics class

    def __init__(self):
    #reads in all configuration files
    #and sets up all mappings

        self.freq_map={'mon':'1mo', 'day':'1d'}


        self.nc_found=[] #list of UM diagnostics found in the NC file during --check_output
        self.nc_missing=[] #list of UM diagnostics missing in the NC file during --check_output

        self.missing=[] # list of diagnostics we failed to add!
        self.added=[] # list of diagnostics we succesfully to added!

        self.nemo_diagnostic_request=[]
        self.nemo_diagnostic_request_off=[]
        self.nemo_diagnostic_request_filename=''
        self.read_ocean_xml()
        self.file_element_id_list=self.get_file_ids()
        #self.rose={}
        
        nemo_field_def_file=main_config['main']['nemo_def']
        if not os.path.isfile(nemo_field_def_file):
            print(nemo_field_def_file+" does not exist")
            exit()
        self.nemo_full_diagnostics = ET.parse(nemo_field_def_file)
        


   
    def get_file_ids(self):
        #compile list of the file element ids
        id_list=[]
        root=self.nemo_diagnostic_request.getroot()
        files=root.findall(".//file")
        if len(files)==0:
            print("XML has no file elements!")
            import pdb; pdb.set_trace()
        for file in files:
            id_list.append(file.attrib['id'])

        return(id_list)
 

        
    def read_rose_app_conf(self,file):
        #read in NEMO configuration XML
        config = configparser.ConfigParser() 
        #turn of the lowercaseization of keys
        config.optionxform = lambda option: option
        configFilePath = file
        if not os.path.isfile(file):
            print("ROSE conf file "+file+" does not exist")
            exit()
           
        with open(configFilePath) as stream:
            config_string=stream.read()

            config.read_string(config_string)
            #config.read_string("[DUMMY]\n" + stream.read())
  
        return(config)

    def read_ocean_xml(self):
        um_nemo_conf="app/xml/rose-app.conf"

        self.rose_conf_file=main_config['user']['job_path']+um_nemo_conf
        self.rose=self.read_rose_app_conf(self.rose_conf_file)

        #find keys containing the ocean diag filename
        diag_keys=[key for key in self.rose.keys() if 'iodef_nemo.xml' in key]
        if len(diag_keys)>1:
            #too many matches!
            print("too many matches?")
            print(diag_keys)
            import pdb; pdb.set_trace()
        #extract link to actual diag file
        nemo_diagnostic_request_file=self.rose[diag_keys[0]]['source'].split('/')[-1]
        this_dir='/'.join(self.rose_conf_file.split('/')[0:-1])
        self.ocean_xml_filename=this_dir+'/file/'+nemo_diagnostic_request_file
        if not os.path.exists(self.ocean_xml_filename):
            print(self.ocean_xml_filename+' does not exist?')
            import pdb; pdb.set_trace()
        #read xml file
    
        self.nemo_diagnostic_request=ET.ElementTree(file=self.ocean_xml_filename)
        #root=tree.getroot()
        self.nemo_diagnostic_request_off=self.get_nemo_commented_fields(self.nemo_diagnostic_request)
        self.nemo_diagnostic_request_filename='app'+self.ocean_xml_filename.split('app')[-1]
        return()



    def get_nemo_commented_fields(self,xml):
        '''
        extract all the commented out diagnostics from a iodef_nemo.xml file
        '''

        sections=[]
        comments=xml.xpath('//comment()')
        #extract all the comments containing fields
        fields=[key.text for key in comments if 'field' in key.text]
        for field in fields:
            sections.append(ET.fromstring(field))
        return(sections)



    
    def nc_check_ocean(self,diag,freq,dims):
        print("Check Ocean output")

        if diag=='ficeberg' and 'olevel' in dims:
            print("ficeberg is actually a 2d field - adjusting")
            dims=dims.replace('olevel','')
        
        if diag=='vowflisf' and not 'olevel' in dims:
            print("vowflisf is actually written out on ocean levels - adjusting")
            dims=dims+' olevel'
        matches=[key for key in nc_output if key.nc_get_variable()==diag]
        if matches:
            print(diag+" found in NC output")
            #import pdb; pdb.set_trace()
            #does this have the required domain?

            #this is slightly odd - but basin is not a dimension, and these basin indices retain longitude for some reason!
            spatial_domain_cf_adjusted=dims.replace('basin','longitude')
            spatial_domain_cf_list=sorted(spatial_domain_cf_adjusted.split())

            #spatial_domain_cf_list=sorted(dims.split())
            ##HERE
            ## SPATIAL domain doesn't map perfectly as ICE is IJ not long lat!
            for this_match in matches:
                nc_domain1=[item.nc_get_variable() for item in this_match.coords().values()]
                nc_domain=sorted([item.replace('nav_lon','longitude').replace('nav_lat','latitude').replace('time_counter','time').replace('deptht','olevel').replace('depthu','olevel').replace('depthv','olevel').replace('depthw','olevel')  for item in nc_domain1])
                #nc_domain=sorted([item.standard_name for item in this_match.coords().values()])
                if spatial_domain_cf_list==nc_domain:
                    #spatial domains match!
                    ##THIS DOESN'T work for CICE
                    #try and get the frequency of this variable from the original file name - should be 1d or 1m
                    match_freq=this_match.properties()['name'].split('_')[1]
                    if freq == match_freq:
                        #time domains match
                        print("Time and spatial domains match")
                        self.nc_found.append([diag,freq,dims])
                        return(True)

        else:
            print(diag+" not found in NC output")

        print("Couldn't find "+diag+" in NC for "+dims+" at "+freq)
        self.nc_missing.append([diag,freq,dims])
        return(False)


    def addOceanDiag(self,diag,freq,dims):
        '''
        add an ocean diagnostic to the ocean XML
        '''

        if check_output:
            self.nc_check_ocean(diag,freq,dims)
            return()
                
        this_freq=self.freq_map[freq]




        this_name_suffix=self.get_name_suffix(diag)
        if this_name_suffix==None:

            print("Couldn't find the name_suffix!")
            print("Unable to add "+diag)
            self.missing.append(diag)
            #import pdb; pdb.set_trace()
            return()

        root=self.nemo_diagnostic_request.getroot()
        fields=root.findall(".//field[@name='"+diag+"']")
        if len(fields)==0:
            plog(diag+" not found in XML diags?")
            #need to pull out of field_def.xml
            new_field=self.get_diag_from_field_def(diag)
            #find the file_group at the correct frequency
            file_group=self.get_file_group(this_freq)
            self.add_field_to_file_group(new_field,file_group,this_name_suffix)
            return()
        else:
            #at least one instance of diagnostic exists in
            #loop over these elements
            for field in fields:
                #get the parent <file> element
                parent=field.getparent()
                #does the parent output_freq equal our required frequency?
                if not 'output_freq' in parent.attrib:
                   print("!!")
                   import pdb; pdb.set_trace()
                   
                if parent.attrib['output_freq']==this_freq:
                    #diagnostic already exists at the requested output frequency
                    #and presumably at the correct grid/name_suffix
                    plog(diag+" already exists at "+freq)
                    return()

            #we are still here, so diag exists, but not at the correct frequency
            plog(diag+" exists, but not at the requested output frequency "+freq)
            #the parent has the correct name_suffix for this diagnostic
            this_name_suffix=parent.attrib['name_suffix']
            #does a <file_group> element exist at the required frequency?
            file_group=self.get_file_group(this_freq)
            #file_groups=root.findall(".//file_group[@output_freq='"+this_freq+"']")
            #if len(file_groups)==0:
            #    #no file group exists for this output freq! Not sure what to do now!
            #    print("No file_group element exists in XML for output frequency "+freq)
            #    import pdb; pdb.set_trace()
            #if len(file_groups)>1:
            #    print("Found more than one file_group for "+freq+" ???")
            #    import pdb; pdb.set_trace()
            #file_group=file_groups[0]

            #file_group,this_name_suffix,this_freq,parent,freq,diag
            #get the <file> elements in this file_group
            files=file_group.findall(".//file")
            if len(files)==0:
                #ah - there are no <file> elements in this file_group!
                plog("No file elements in this file_group!")
                #need to copy one across
                #create new file element in this file_group
                new_file_element=ET.SubElement(file_group,'file')
                #copy attributes and text  from field parent (file elements)
                #create a dictionary from the parent attributes
                parent_attrib=dict(parent.attrib)
                #change to correct output_freq
                parent_attrib['output_freq']=this_freq
                #modify filename prefix
                parent_attrib['name']='@expname@_'+freq
                #create new file id
                file_id=self.get_unique_file_id()
                plog("Creating new file element "+file_id+" for "+freq)
                parent_attrib['id']=file_id
                self.file_element_id_list.append(file_id)
                #set the updated attributes#
                new_file_element.attrib.update(parent_attrib)
                #set the text (if there is any)
                new_file_element.text=parent.text
                plog("Adding "+diag+" to "+file_id)
                new_file_element.append(deepcopy(field))
                self.added.append(diag)
                print("Done")
                return()

            #here - there are multiple files in this fille_group for the chosen output freq
            #we need to check to see if one matches the grid (or name_suffix) of the requested diagnostic
            #name_suffix_found=False
            for file in files:
                if file.attrib['name_suffix']==this_name_suffix:
                    #found the grid
                    #this file is the correct place to add the diagnostic
                    plog("copying "+diag+" to "+file.attrib['id'])
                    self.added.append(diag)
                    file.append(deepcopy(fields[0]))
                    return()

            print("Couldn't find "+this_name_suffix+"in "+file_groups[0])
            import pdb; pdb.set_trace()

       
     
    def addOceanDiag_old2(self,line):

        diag=line['variable']
        freq=line['time']
        dims=line['space']

        #print(+" "+freq)

        #diag = cf variable name
        #freq= output frequency mon,day etc
        print("")

        print("-------------------------")
        print("Ocean: "+diag)

        if diag=='areacello':
           #this is just the ocean grid cell area - not a model diagnostic
           plog(diag+" is just the ocean grid cell area - don't need to write this out as a diagnostics")
           return()

        #is this really a CICE diagnostic?

        cice_diagnostics=cice.rose['namelist:icefields_nml']
        #if this diag is not mapped in the cf mapping AND exists in the cice_diagnostics
        #it must be a cice diagnostic!
        if not diag in um.cf_to_stash and 'f_'+diag in cice_diagnostics:
           #YES!
           print(diag+" is a CICE diagnostic - adding to CICE")
           cice.addIceDiag(line)
           return()
        
        
        #is this diag mapped in the cf mapping table?
        if diag in um.cf_to_stash:
           cf_mapping_expression=um.cf_to_stash[diag]['expression']
           if not cf_mapping_expression==diag:
              #the mapping expression is not just the diag name
              #need to break down the mapping expression and add each sub diagnostic in turn
              #these could be atmosphere(UM) diags or NEMO diags or even CICE diags!
              plog(diag+" is mapped in the cf mappings")
              print(cf_mapping_expression)

              self.add_diags_in_mapping(cf_mapping_expression,line)
              return()
           #here - the mapping expresion is JUST the diag (this can happen if a mask is being applied etc)
           
           
                 
                 
              # we need to extract JUST the diagnostics
              #pattern for removing functions of the form aaaa_bbbb(
              pattern=r'[0-9a-z_]+\('
              #remove pattern and any \n
              cf_mapping_expression2=re.sub(pattern,'',cf_mapping_expression.replace('\n',' '))
              #translation_table = str.maketrans("", "", ",*+-/().0123456789_ABCDEFGHIJKLMNOPQRSTUVWXYZ")
              translation_table = str.maketrans("", "", ",+/()")
              #strip out all operators, numbers and upper case characters
              sub_diag_string=cf_mapping_expression2.translate(translation_table)
              #split this into a list of nemo diagnostics, excluding empty strings
              nemo_diags=[x for x in sub_diag_string.split(' ') if x !='']
              print("------------")
              pattern = re.compile('[^a-z_]+')

              nemo_diags_all=[]
              for nemo_diag in nemo_diags:
                 #exclude upper case constants
                 #exclude names including "_0" these are reference diagnostics for calc_zostoga() - but come from the PI, not this job
                 # see https://code.metoffice.gov.uk/trac/cdds/browser/main/trunk/mip_convert/mip_convert/process/processors.py
                 if not (any(char.isupper() for char in nemo_diag) or nemo_diag=='*' or '=' in nemo_diag or nemo_diag=='-1' or nemo_diag=='-' or '_0' in nemo_diag):
                    n1=re.sub(pattern,'',nemo_diag)
                    nemo_diags_all.append(n1)

              nemo_diags_unique=list(set(nemo_diags_all))

              diag_in_mapping=False
              if diag in nemo_diags_unique:
                 #the diag reappears in the mapping!
                 # this can happen for situations like:
                 #correct_evaporation(evs, sowaflup - (evs - (ficeberg+ friver + prsn + pr + fsitherm) ), areacello)
                 #in this case we will exit to the remainder of the function once the other diagnostics in the mapping have been dealt with
                 #this will then add the original diagnostics diag
                 diag_in_mapping=True
              
              if len(nemo_diags_unique)==1 and nemo_diags_unique[0]==diag:
                 #The mapping file just points back to the diag!
                 #something lime mask_copy(tauuo, mask_2D_U)
                 #we don't care about the masking here - just need to output the correct diagnostics
                 print("Mapping just points back to this diagnostics, so don't need to worry about this")
                 import pdb; pdb.set_trace()
                 #UMO why doesn't that work?
              else:
                 ###HERE
                 ## need to filter things like so_0 ??
                 #Then loop over all and call addOceanDiag again
                 plog(diag+" is mapped to the following diagnostics in the mapping conf: "+' '.join(nemo_diags_unique))
                 plog("Adding these diagnostics..")
                 

                 for this_diag in nemo_diags_unique:
                    if this_diag!=diag and this_diag!='':
                       #only proceed if this is not the same as the diag requested!
                       #to avoid infinite loop
                       #copy the line sent to addOceanDiag()
                       this_line=dict(line)
                       this_line['variable']=this_diag
                    
                       self.addOceanDiag(this_line)
                 if not diag_in_mapping:      
                    print("All added")
                    return()
                 #we still need to add the original diag - so carry on a do that
                 #import pdb; pdb.set_trace()


        else:         
           #diag does NOT appear in the current um.cf_to_stash      

           #freq_map={'mon':'1mo', 'day':'1d'}
           this_freq=self.freq_map[freq]
           this_name_suffix=self.get_name_suffix(diag)
           if this_name_suffix==None:

               print("Couldn't find the name_suffix!")
               print("Unable to add "+diag)
               self.missing.append(diag)
               #import pdb; pdb.set_trace()
               return()

           root=self.nemo_diagnostic_request.getroot()
           fields=root.findall(".//field[@name='"+diag+"']")
           if len(fields)==0:
               plog(diag+" not found in XML diags?")
               #need to pull out of field_def.xml
               new_field=self.get_diag_from_field_def(diag)
               #find the file_group at the correct frequency
               file_group=self.get_file_group(this_freq)
               self.add_field_to_file_group(new_field,file_group,this_name_suffix)
               return()
           else:
               #at least one instance of diagnostic exists in
               #loop over these elements
               for field in fields:
                   #get the parent <file> element
                   parent=field.getparent()
                   #does the parent output_freq equal our required frequency?
                   if parent.attrib['output_freq']==this_freq:
                       #diagnostic already exists at the requested output frequency
                       #and presumably at the correct grid/name_suffix
                       plog(diag+" already exists at "+freq)
                       return()

               #we are still here, so diag exists, but not at the correct frequency
               plog(diag+" exists, but not at the requested output frequency "+freq)
               #the parent has the correct name_suffix for this diagnostic
               this_name_suffix=parent.attrib['name_suffix']
               #does a <file_group> element exist at the required frequency?
               file_group=self.get_file_group(this_freq)
               #file_groups=root.findall(".//file_group[@output_freq='"+this_freq+"']")
               #if len(file_groups)==0:
               #    #no file group exists for this output freq! Not sure what to do now!
               #    print("No file_group element exists in XML for output frequency "+freq)
               #    import pdb; pdb.set_trace()
               #if len(file_groups)>1:
               #    print("Found more than one file_group for "+freq+" ???")
               #    import pdb; pdb.set_trace()
               #file_group=file_groups[0]

               #file_group,this_name_suffix,this_freq,parent,freq,diag
               #get the <file> elements in this file_group
               files=file_group.findall(".//file")
               if len(files)==0:
                   #ah - there are no <file> elements in this file_group!
                   plog("No file elements in this file_group!")
                   #need to copy one across
                   #create new file element in this file_group
                   new_file_element=ET.SubElement(file_group,'file')
                   #copy attributes and text  from field parent (file elements)
                   #create a dictionary from the parent attributes
                   parent_attrib=dict(parent.attrib)
                   #change to correct output_freq
                   parent_attrib['output_freq']=this_freq
                   #modify filename prefix
                   parent_attrib['name']='@expname@_'+freq
                   #create new file id
                   file_id=self.get_unique_file_id()
                   plog("Creating new file element "+file_id+" for "+freq)
                   parent_attrib['id']=file_id
                   self.file_element_id_list.append(file_id)
                   #set the updated attributes#
                   new_file_element.attrib.update(parent_attrib)
                   #set the text (if there is any)
                   new_file_element.text=parent.text
                   plog("Adding "+diag+" to "+file_id)
                   new_file_element.append(deepcopy(field))
                   self.added.append(diag)
                   print("Done")
                   return()

               #here - there are multiple files in this fille_group for the chosen output freq
               #we need to check to see if one matches the grid (or name_suffix) of the requested diagnostic
               #name_suffix_found=False
               for file in files:
                   if file.attrib['name_suffix']==this_name_suffix:
                       #found the grid
                       #this file is the correct place to add the diagnostic
                       plog("copying "+diag+" to "+file.attrib['id'])
                       self.added.append(diag)
                       file.append(deepcopy(fields[0]))
                       return()

               print("Couldn't find "+this_name_suffix+"in "+file_groups[0])
               import pdb; pdb.set_trace()



    def add_diags_in_mapping(self,expression,line):
       '''
       breaks down the expression into UM, NEMO and CICE diagnostics
       and adds them to the relevant diagnostic list using the freq and space information in line
       '''
       freq=line['time']
       dims=line['space']

       ##UM diagnostics
       #does the expression contain any UM stash diags?
       pattern_s = r'm\d{2}s\d{2}i\d{3}'
       matches = re.findall(pattern_s, cf_mapping_expression)
       if matches:
          #this is a stash diagnostic
          plog(cf_mapping_expression+" contains UM stash diagnostics")
          plog("adding stash diagnostics")
          for match in matches:
             new_line=dict(line)
             #new_line['variable']=
             um.add_stash(match,freq,dims)
             import pdb; pdb.set_trace()
       
           


    def addOceanDiag_old(self,line):

        diag=line['variable']
        freq=line['time']
        dims=line['space']

        #print(+" "+freq)

        #diag = cf variable name
        #freq= output frequency mon,day etc
        print("")

        print("-------------------------")
        print("Ocean: "+diag)

        if diag=='areacello':
           #this is just the ocean grid cell area - not a model diagnostic
           plog(diag+" is just the ocean grid cell area - don't need to write this out as a diagnostics")
           return()

        #is this really a CICE diagnostic?

        cice_diagnostics=cice.rose['namelist:icefields_nml']
        #if this diag is not mapped in the cf mapping AND exists in the cice_diagnostics
        #it must be a cice diagnostic!
        if not diag in um.cf_to_stash and 'f_'+diag in cice_diagnostics:
           #YES!
           print(diag+" is a CICE diagnostic - adding to CICE")
           cice.addIceDiag(line)
           return()
        
        
        #is this diag mapped in the cf mapping table?
        if diag in um.cf_to_stash:
           cf_mapping_expression=um.cf_to_stash[diag]['expression']
           if not cf_mapping_expression==diag:
              #the mapping expression is not just the diag name
              plog(diag+" is mapped in the cf mappings")
              print(cf_mapping_expression)
              #check to see if this is really a UM stash diagnostic!
              pattern_s = r'm\d{2}s\d{2}i\d{3}\[.*?\]'
              matches = re.findall(pattern_s, cf_mapping_expression)
              if matches:
                 #this is a stash diagnostic
                 plog(cf_mapping_expression+" is a stash diagnistics for the UM")
                 plog("adding stash diagnostics")
                 um.add_cf_diagnostic(line)
                 return()
                 
                 
              # we need to extract JUST the diagnostics
              #pattern for removing functions of the form aaaa_bbbb(
              pattern=r'[0-9a-z_]+\('
              #remove pattern and any \n
              cf_mapping_expression2=re.sub(pattern,'',cf_mapping_expression.replace('\n',' '))
              #translation_table = str.maketrans("", "", ",*+-/().0123456789_ABCDEFGHIJKLMNOPQRSTUVWXYZ")
              translation_table = str.maketrans("", "", ",+/()")
              #strip out all operators, numbers and upper case characters
              sub_diag_string=cf_mapping_expression2.translate(translation_table)
              #split this into a list of nemo diagnostics, excluding empty strings
              nemo_diags=[x for x in sub_diag_string.split(' ') if x !='']
              print("------------")
              pattern = re.compile('[^a-z_]+')

              nemo_diags_all=[]
              for nemo_diag in nemo_diags:
                 #exclude upper case constants
                 #exclude names including "_0" these are reference diagnostics for calc_zostoga() - but come from the PI, not this job
                 # see https://code.metoffice.gov.uk/trac/cdds/browser/main/trunk/mip_convert/mip_convert/process/processors.py
                 if not (any(char.isupper() for char in nemo_diag) or nemo_diag=='*' or '=' in nemo_diag or nemo_diag=='-1' or nemo_diag=='-' or '_0' in nemo_diag):
                    n1=re.sub(pattern,'',nemo_diag)
                    nemo_diags_all.append(n1)

              nemo_diags_unique=list(set(nemo_diags_all))

              diag_in_mapping=False
              if diag in nemo_diags_unique:
                 #the diag reappears in the mapping!
                 # this can happen for situations like:
                 #correct_evaporation(evs, sowaflup - (evs - (ficeberg+ friver + prsn + pr + fsitherm) ), areacello)
                 #in this case we will exit to the remainder of the function once the other diagnostics in the mapping have been dealt with
                 #this will then add the original diagnostics diag
                 diag_in_mapping=True
              
              if len(nemo_diags_unique)==1 and nemo_diags_unique[0]==diag:
                 #The mapping file just points back to the diag!
                 #something lime mask_copy(tauuo, mask_2D_U)
                 #we don't care about the masking here - just need to output the correct diagnostics
                 print("Mapping just points back to this diagnostics, so don't need to worry about this")
                 import pdb; pdb.set_trace()
                 #UMO why doesn't that work?
              else:
                 ###HERE
                 ## need to filter things like so_0 ??
                 #Then loop over all and call addOceanDiag again
                 plog(diag+" is mapped to the following diagnostics in the mapping conf: "+' '.join(nemo_diags_unique))
                 plog("Adding these diagnostics..")
                 

                 for this_diag in nemo_diags_unique:
                    if this_diag!=diag and this_diag!='':
                       #only proceed if this is not the same as the diag requested!
                       #to avoid infinite loop
                       #copy the line sent to addOceanDiag()
                       this_line=dict(line)
                       this_line['variable']=this_diag
                    
                       self.addOceanDiag(this_line)
                 if not diag_in_mapping:      
                    print("All added")
                    return()
                 #we still need to add the original diag - so carry on a do that
                 #import pdb; pdb.set_trace()
                 
        #diag does NOT appear in the current um.cf_to_stash      

        #freq_map={'mon':'1mo', 'day':'1d'}
        this_freq=self.freq_map[freq]
        this_name_suffix=self.get_name_suffix(diag)
        if this_name_suffix==None:
            
            print("Couldn't find the name_suffix!")
            print("Unable to add "+diag)
            self.missing.append(diag)
            #import pdb; pdb.set_trace()
            return()

        root=self.nemo_diagnostic_request.getroot()
        fields=root.findall(".//field[@name='"+diag+"']")
        if len(fields)==0:
            plog(diag+" not found in XML diags?")
            #need to pull out of field_def.xml
            new_field=self.get_diag_from_field_def(diag)
            #find the file_group at the correct frequency
            file_group=self.get_file_group(this_freq)
            self.add_field_to_file_group(new_field,file_group,this_name_suffix)
            return()
        else:
            #at least one instance of diagnostic exists in
            #loop over these elements
            for field in fields:
                #get the parent <file> element
                parent=field.getparent()
                #does the parent output_freq equal our required frequency?
                if parent.attrib['output_freq']==this_freq:
                    #diagnostic already exists at the requested output frequency
                    #and presumably at the correct grid/name_suffix
                    plog(diag+" already exists at "+freq)
                    return()

            #we are still here, so diag exists, but not at the correct frequency
            plog(diag+" exists, but not at the requested output frequency "+freq)
            #the parent has the correct name_suffix for this diagnostic
            this_name_suffix=parent.attrib['name_suffix']
            #does a <file_group> element exist at the required frequency?
            file_group=self.get_file_group(this_freq)
            #file_groups=root.findall(".//file_group[@output_freq='"+this_freq+"']")
            #if len(file_groups)==0:
            #    #no file group exists for this output freq! Not sure what to do now!
            #    print("No file_group element exists in XML for output frequency "+freq)
            #    import pdb; pdb.set_trace()
            #if len(file_groups)>1:
            #    print("Found more than one file_group for "+freq+" ???")
            #    import pdb; pdb.set_trace()
            #file_group=file_groups[0]

            #file_group,this_name_suffix,this_freq,parent,freq,diag
            #get the <file> elements in this file_group
            files=file_group.findall(".//file")
            if len(files)==0:
                #ah - there are no <file> elements in this file_group!
                plog("No file elements in this file_group!")
                #need to copy one across
                #create new file element in this file_group
                new_file_element=ET.SubElement(file_group,'file')
                #copy attributes and text  from field parent (file elements)
                #create a dictionary from the parent attributes
                parent_attrib=dict(parent.attrib)
                #change to correct output_freq
                parent_attrib['output_freq']=this_freq
                #modify filename prefix
                parent_attrib['name']='@expname@_'+freq
                #create new file id
                file_id=self.get_unique_file_id()
                plog("Creating new file element "+file_id+" for "+freq)
                parent_attrib['id']=file_id
                self.file_element_id_list.append(file_id)
                #set the updated attributes#
                new_file_element.attrib.update(parent_attrib)
                #set the text (if there is any)
                new_file_element.text=parent.text
                plog("Adding "+diag+" to "+file_id)
                new_file_element.append(deepcopy(field))
                self.added.append(diag)
                print("Done")
                return()

            #here - there are multiple files in this fille_group for the chosen output freq
            #we need to check to see if one matches the grid (or name_suffix) of the requested diagnostic
            #name_suffix_found=False
            for file in files:
                if file.attrib['name_suffix']==this_name_suffix:
                    #found the grid
                    #this file is the correct place to add the diagnostic
                    plog("copying "+diag+" to "+file.attrib['id'])
                    self.added.append(diag)
                    file.append(deepcopy(fields[0]))
                    return()

            print("Couldn't find "+this_name_suffix+"in "+file_groups[0])
            import pdb; pdb.set_trace()

    def get_name_suffix(self,diag):
        #returns the correct name_suffix/grid for a given diagnostic
        #diag is a cf variable name
        #need to work back to a field_ref and then cross ref with a field id to find the field group parent
        #this parent id should be the name suffix
        root=self.nemo_diagnostic_request.getroot()
        fields=root.findall(".//field[@name='"+diag+"']")
        if len(fields)==0:
            #we didn't find anything!
            #try field_def.xml
            fields=self.nemo_full_diagnostics.findall(".//field[@name='"+diag+"']")
            if len(fields)==0:
                #nothing here either!
                plog(bold(diag+" is unknown cf variable!"))
                #loop over extra nemo fields added in main_config
                for i in [key for key in main_config if 'nemo_field' in key]:
                   if diag in i:
                      print("Mapping for "+diag+" found in config!")
                      this_nemo_field_ref=main_config[i]['field_ref']
                      this_nemo_file=main_config[i]['nemo_file']
                      print("Adding "+this_nemo_field_ref+" to the nemo "+this_nemo_file+" file stream")
                      ##Can we add this directly? But really need to make new mapping..?
                      
                      import pdb; pdb.set_trace()
                      

                #does this
                print("You might have to add this by hand to "+self.nemo_diagnostic_request_filename)
                self.missing.append(diag)
                return(None)

            #Found at least one field with matching name
            #extract the name from the first field attrinb
            field_name=fields[0].attrib['name']
            field_ref=fields[0].attrib['field_ref']
            #find the diagnostics with this field id
            diags=self.nemo_full_diagnostics.findall(".//field[@id='"+field_ref+"']")
            if len(diags)==0:
                print("No diags found with field_ref "+field_ref)
                import pdb; pdb.set_trace()
            if len(diags)>1:
                print("Multiple matching ids!")
                import pdb; pdb.set_trace()
            this_parent=diags[0].getparent()
            field_group_id=this_parent.attrib['id']
            #does this_parent have a grid_ref attrib?
            if not 'grid_ref' in this_parent.attrib:
                #no grid_ref, just return the id as the nam_suffix
                return('_'+field_group_id)
            grid_ref=this_parent.attrib['grid_ref']
            #grid ref should either be "grid_X_YD" X=T,U,V,W Y=2,3
            # OR "scalar"
            if grid_ref=='scalar':
                return('_scalar')
            #must be grid_X_YD?
            if not '_' in grid_ref:
                print("Not sure what grid we have here!")
                print(grid_ref)
                import pdb; pdb.set_trace()
            #return truncated grid ref (removing _2D or _3D)
            return('_'+grid_ref[:-3])

            #import pdb; pdb.set_trace()
            ##### THIS NEEDS FIXING HERE
            ##01/Dec/2023
            #this_diag_attrib=dict(diags[0].attrib)
            #this_diag_attrib['name']=field_name
            #import pdb; pdb.set_trace()



        #find the parent
        parent=fields[0].getparent()
        if parent==None:
            print("No Parent!!")
            import pdb; pdb.set_trace()
        #we will assume the name_suffix is what we need here! '_grid_T' '_grid_U' '_grid_V' '_diaptr'
        if 'name_suffix' in parent.attrib:
            grid=parent.attrib['name_suffix']
            return(grid)
        else:
            print('No name_suffix attribute in the parent!')
            import pdb; pdb.set_trace()


    def get_diag_from_field_def(self,diag):
        fields=self.nemo_full_diagnostics.findall(".//field[@name='"+diag+"']")

        if len(fields)==0:
            #nothing here either!
            plog(diag+" is unknown cf variable!")
            return(None)
        #Found at least one field with matching name
        #extract the name from the first field attrinb
        field_name=fields[0].attrib['name']
        field_ref=fields[0].attrib['field_ref']

        #find the diagnostics with this field id
        diags=self.nemo_full_diagnostics.findall(".//field[@id='"+field_ref+"']")
        if len(diags)==0:
            print("No diags found with field_ref "+field_ref)
            import pdb; pdb.set_trace()
        if len(diags)>1:
            print("Multiple matching ids!")
            import pdb; pdb.set_trace()
        new_diag=deepcopy(diags[0])
        new_diag.attrib['name']=field_name
        #convert the id to field_ref
        new_diag.attrib['field_ref']=new_diag.attrib['id']
        #delete existing id entry
        del new_diag.attrib['id']
        new_cell_measure=deepcopy(self.nemo_diagnostic_request.findall(".//variable[@name='cell_measures']")[0])
        new_diag.append(new_cell_measure)
        return(new_diag)

            

    def add_field_to_file_group(self,field,file_group,name_suffix):
        '''
        
        '''

        #get the <file> elements in this file_group
        root=self.nemo_diagnostic_request.getroot()
        files=file_group.findall(".//file")
        if len(files)==0:
                #ah - there are no <file> elements in this file_group!
                plog("No file elements in this file_group!")
                #import pdb; pdb.set_trace()
                #need to copy one across
                #create new file element in this file_group
                new_file_element=ET.SubElement(file_group,'file')
                #copy attributes and text  from file element that matches name_suffix
                match_files=root.findall(".//file[@name_suffix='"+name_suffix+"']")
                if len(match_files)==0:
                    print("No file elemnts found with a matching name suffix! ")
                    import pdb; pdb.set_trace()
                new_attrib=dict(match_files[0].attrib)

                #change to correct output_freq
                this_freq=file_group.attrib['output_freq']
                freq=None
                for key in self.freq_map:
                    if self.freq_map[key]==this_freq:
                        freq=key
                if freq==None:
                    print("Freq mapping not found for "+this_freq)
                    import pdb; pdb.set_trace()

                #freq=
                new_attrib['output_freq']=file_group.attrib['output_freq']
                #modify filename prefix
                new_attrib['name']='@expname@_'+freq
                #create new file id
                file_id=self.get_unique_file_id()
                plog("Creating new file element "+file_id)
                new_attrib['id']=file_id
                file_element_id_list.append(file_id)
                #set the updated attributes#
                new_file_element.attrib.update(new_attrib)
                #create a dictionary from the parent attributes
                diag=field.attrib['name']
                new_file_element.text="\n"
                plog("Adding "+diag+" to "+file_id)
                new_file_element.append(deepcopy(field))
                self.added.append(diag)
                print("Done")
                return()
        #here - there are multiple files in this fille_group for the chosen output freq
        #we need to check to see if one matches the grid (or name_suffix) of the requested diagnostic
        name_suffix_found=False
        for file in files:
            if file.attrib['name_suffix']==name_suffix:
                #found the grid
                grid_found=True
                diag=field.attrib['name']
                #this file is the correct place to add the diagnostic
                plog("copying "+diag+" to "+file.attrib['id'])
                file.append(deepcopy(field))
                self.added.append(diag)
                return()

        
        print("Couldn't find "+name_suffix+" in "+file_group.attrib['id'])
        
        import pdb; pdb.set_trace()


    def get_file_group(self,this_freq):
        root=self.nemo_diagnostic_request.getroot()
        file_groups=root.findall(".//file_group[@output_freq='"+this_freq+"']")
        if len(file_groups)==0:
            #no file group exists for this output freq! Not sure what to do now!
            print("No file_group element exists in XML for output frequency "+this_freq)
            import pdb; pdb.set_trace()
        if len(file_groups)>1:
            print("Found more than one file_group for "+this_freq+" ???")
            import pdb; pdb.set_trace()
        file_group=file_groups[0]
        return(file_group)

    
    def get_unique_file_id(self):
        #returns a unique id in the format of file_element_id_list, but unique
        #assumes format is fileNNNN where NNNN is a number
        n_list=[]
        for i in self.file_element_id_list:
            n_list.append(int(i[4:]))
        n_list.sort()
        new_index=n_list[-1]+1
        new_file_id='file'+str(new_index)
        return(new_file_id)


    def write(self,output_file):
        self.nemo_diagnostic_request.write(output_file)
        print("Written "+output_file)

    
    pass



    
#freq_map={'mon':'1mo', 'day':'1d'}
#
#space_mappings={'longitude latitude time':"'DIAG'",
#                'longitude latitude height2m time':"'DIAG'",
#                'longitude latitude height10m time':"'DIAG'",
#                'longitude latitude alevel time':"'DALLTH'",
#                'longitude latitude alevhalf time':"'DALLRH'",
#                'longitude latitude plev8 time':"'PLEV8'",
#                'longitude latitude plev19 time':"'PLEV19'"
#                 }
#
#freq_mappings={'mon':"'TMONMN'",
#               'day':"'TDAYMN'"
#               }
#
#
#default_usage={'mon':"'UPM'",
#               'day':"'UPD'"
#               }
#



def read_cf_diagnostics():
    cf_diagnostics_file=main_config['user']['cf_diagnostics_file']

    if not os.path.isfile(cf_diagnostics_file):
        print("CF Diagnostics CSV file "+cf_diagnostics_file+" does not exist")
        exit()

    # Open the CSV cf_diagnostics file
    with open(cf_diagnostics_file, 'r') as infile:
        # Create a CSV reader object for reading as a dictionary
        csv_reader = csv.DictReader(infile)

        # Convert the CSV data into a list of dictionaries
        variable_list = list(csv_reader)
    return(variable_list)


def read_config(conf_file):
    if not os.path.isfile(conf_file):
        print("Config file "+conf_file+" dos not exist")
        exit()

    main_config = configparser.ConfigParser()   
    #turn of the lowercaseization of keys
    main_config.optionxform = lambda option: option
    main_config.read(conf_file)
    return(main_config)

def start_logging():
    if 'log_file' in main_config['user']:
        log_file=main_config['user']['log_file'].strip("'")
    else:
        log_file='cf_to_um.log'    
    #logging.basicConfig(filename=log_file, encoding='utf-8', level=logging.INFO, filemode='w')
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode='w')
    print("Writing log to "+log_file)

def check_stash_eq(um,nemo):
    #checks that the stash entries in um.rose and nemo.rose are identical!¬
    #They SHOULD be!
    nemo_keys=[key for key in nemo.rose.keys() if 'umstash' in key]
    um_keys=[key for key in um.rose.keys() if 'umstash' in key]
    import pdb; pdb.set_trace()


def plog(message):
    print(message)
    logging.info(message)
    
    

##################################

parser = argparse.ArgumentParser(description='cf_to_um_diagnostics.py adds CF diagnostics to existing UM rose job')
parser.add_argument('-c', '--config',type=str)   
parser.add_argument('-s', '--stash',type=str,choices=['um','xios'])
parser.add_argument('-z', '--check_output')

args = parser.parse_args()

if args.config:
    conf_file=args.config
else:
    conf_file='cf_to_um.conf'

if args.stash:
    stash_type=args.stash
else:
    stash_type='um'

#this option allows use to check the netcdf/pp output 
check_output=False
if args.check_output:
    print("Checking NC output..")
    import cf
    nc_output=cf.read(args.check_output+"/*nc")
    check_output=True

#read in main config file
main_config=read_config(conf_file)


#setup logging file
start_logging()
#Now rose should have all the required time and space domains defined as in the freq_mappings and space_mappings
plog("Adding to the "+bold(stash_type)+" STASH, Nemo and CICE diagnostics")
print("------------")

#read in cf mappings
cf_mappings=read_cf_mappings()

#read in cf variable list
variable_list=read_cf_diagnostics()



#initialize um stash/cice instance
#stash_type is which STASH to add diagnostics to UM or XML
um=UM(stash_type)

#initialize Nemo instance
nemo=Nemo()
#initialize CICE instance
cice=CICE()

#if not check_stash_eq(um,nemo):
#    print("XML and UM STASH entries differ!")
#    import pdb; pdb.set_trace()

   

#---------------------





print("----------------------------")
#loop over all cf variables

#UM has multiple realms
um_realms=['atmos','landIce','land']

for line in variable_list:
    diag=line['variable']
    freq=line['time']
    dims=line['space']
    add_cf_diagnostic(diag,freq,dims)
    #import pdb; pdb.set_trace()

#    realm=line['realm']
#    #if realm in um_realms:
#    if any(element in um_realms for element in realm.split(' ')):
#        #all realms here use STASH
#        print(".")
#        um.add_cf_diagnostic(line)
#    elif 'ocean' in realm:
#        #add ocean diagnostic
#        nemo.addOceanDiag(line)
#        #print()
#    elif 'seaIce' in realm:
#        #add Sea Ice diagnostic
#        cice.addIceDiag(line)
#        #print()
#    else:
#        print("Unknown Realm?")
#        print(line['realm'])
#        import pdb; pdb.set_trace()

if check_output:
    print("")
    print("")
    if um.nc_found:
        print("The following STASH diagnostics were found in the output")
        for i in um.nc_found:
            print(f'{i[0]:10}  {i[1]} [{i[2]}] ')
        print("--------------------------")
        
    if nemo.nc_found:
        print("The following NEMO diagnostics were found in the output")
        for i in nemo.nc_found:
            print(f'{i[0]:20}  {i[1]} [{i[2]}] ')
        print("--------------------------")

    if cice.nc_found:
        print("The following CICE diagnostics were found in the output")
        for i in cice.nc_found:
            print(f'{i[0]:10}  {i[1]} [{i[2]}] ')
        print("--------------------------")

    print("")
            
    if um.nc_missing:
        print("The following STASH diagnostics are "+color.BOLD+" missing"+color.END+"from the output")
        for i in um.nc_missing:
            print(f'{i[0]:10}  {i[1]} [{i[2]}] ')
        print("--------------------------")
    else:
        print("There were no missing STASH diagnostics")

    if nemo.nc_missing:
        print("The following NEMO diagnostics are "+color.BOLD+" missing"+color.END+" from the output")
        for i in nemo.nc_missing:
            print(f'{i[0]:20}  {i[1]} [{i[2]}] ')
        print("--------------------------")
    else:
        print("There were no missing NEMO diagnostics")

    if cice.nc_missing:
        print("The following CICE diagnostics are "+color.BOLD+" missing"+color.END+" from the output")
        for i in cice.nc_missing:
            print(f'{i[0]:10}  {i[1]} [{i[2]}] ')
        print("--------------------------")
    else:
        print("There were no missing CICE diagnostics")

    exit()        



um.missing.sort()
um.added.sort()
nemo.missing.sort()
nemo.added.sort()
cice.missing.sort()
cice.added.sort()
        
plog(bold("UM diagnostics unable to add: "+' '.join(um.missing)))
plog(bold("Nemo diagnostics unable to add: "+' '.join(nemo.missing)))
plog(bold("CICE diagnostics unable to add: "+' '.join(cice.missing)))
print("----------------------")
plog(bold("UM diagnostics added: "+' '.join(um.added)))
plog(bold("NEMO diagnostics added: "+' '.join(nemo.added)))
plog(bold("CICE diagnostics added: "+' '.join(cice.added)))

#write diagnostics definition files
#extract job name from filename
pattern = r'/u-([a-zA-Z0-9-]+)/'
jobname=re.search(pattern,um.rose_stash).group(1)
#jobname=re.search(pattern,main_config['user']['job_path']).group(1)
ice_conf="app/nemo_cice/rose-app.conf"
um_nemo_conf="app/xml/rose-app.conf"
#ice_output_filename=jobname+'_'+ice_conf.replace('/','__')
### NOT OUTPUTTING CORRECT UM CONF!!!
um_output_filename=um.rose_stash.split('roses/')[-1].replace('/','__')
ocean_output_filename=nemo.ocean_xml_filename.split('roses/')[-1].replace('/','__')
cice_output_filename=cice.rose_cice.split('roses/')[-1].replace('/','__')

um.write(um_output_filename)
nemo.write(ocean_output_filename)
cice.write(cice_output_filename)
exit()    
           
# LBPROC Processing code. This indicates what processing has been done to the basic field. It should be
#0 if no processing has been done, otherwise the relevant numbers from the list below are added together.
#E.g. for a field of minimum cloudbase during a time period, LBPROC=4096, and the value will be ”missing
#data” for any point where at least one timestep in the period has MDI, indicating intermittent cloud. But
#if the UM STASH domain profile masking imsk=4, LBPROC=4096+252144=256240, and the value at a
#point will only be ”missing data” if there is never any cloud present during the time period.
#1 Difference from another experiment.
#2 Difference from zonal (or other spatial) mean.
#4 Difference from time mean.
#8 X-derivative (d/dx)
#16 Y-derivative (d/dy)
#32 Time derivative (d/dt)
#64 Zonal mean field
#128 Time mean field or accumulation
#256 Product of two fields
#512 Square root of a field
#23 c© Crown Copyright 2023
#UMDP: F03
#Input and Output File Formats
#Chapter 4. LOOKUP table and PP headers
#1024 Difference between fields at levels BLEV and BRLEV
#2048 Mean over layer between levels BLEV and BRLEV
#4096 Minimum value of field during time period
#8192 Maximum value of field during time period
#16384 Magnitude of a vector, not specifically wind speed
#32768 Log10 of a field
#65536 Variance of a field
#131072 Mean over an ensemble of paralle
