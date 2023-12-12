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
import configparser
import argparse
import logging
import uuid
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

def bold(message):
    return(color.BOLD+message+color.END)

 
#############  Atmosphere/Land class for UM
class UM:
    #UM stash  class
    # add_cf_diagnostic() adds an atmosphere/land/landice cf variable as the require STASH codes
    def __init__(self,umOrXIOS):
        #we can use the STASH in the um app or the STASH in the xml app (for netcdf)
        #umOrXIOS = 'um' or 'xios'
        #reads in all configuration files
        #and sets up all mappings

        self.missing=[] # list of diagnostics we failed to add!
        self.added=[] # list of diagnostics we succeeded in adding!
        self.umOrXIOS=umOrXIOS
        self.default_usage={'mon':"'UPM'",
                            'day':"'UPD'"
                            }

        
        self.freq_mappings={'mon':"'TMONMN'",
                            'day':"'TDAYMN'"
                            }
        self.space_mappings={'longitude latitude time':"'DIAG'",
                             'longitude latitude height2m time':"'DIAG'",
                             'longitude latitude height10m time':"'DIAG'",
                             'longitude latitude alevel time':"'DALLTH'",
                             'longitude latitude alevhalf time':"'DALLRH'",
                             'longitude latitude plev8 time':"'PLEV8'",
                             'longitude latitude plev19 time':"'PLEV19'",
                             'longitude latitude soil time':"'DSOIL'",
                             'longitude latitude sdepth time':"'DSOIL'",
                             'longitude latitude sdepth1 time':"'DSOIL1'"
                             }

        self.cf_to_stash={}
        self.rose_time_domain_mappings={}
        self.rose_space_domain_mappings={}
        self.use_matrix={}
        self.stash_levels={}
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
        cmip6_rose=main_config['main']['cmip6']
        #read in standard CMIP6 time and spatial domain definitions
        self.cmip6,self.cmip6_header=self.read_rose_app_conf(cmip6_rose)
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
        
        #Now rose should have all the required time and space domains defined as in the freq_mappings and space_mappings


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
                bits=line.split('|')
                if bits[0]=='1':
                    model=bits[1].strip(' ')
                    if model=='-1':
                        #end of file
                        break
                    sec=bits[2].strip(' ')
                    item=bits[3].strip(' ')
                    name=bits[4]
                    scode='m0'+model+'s'+sec.zfill(2)+'i'+item.zfill(3)
                if bits[0]=='2':
                    level=int(bits[5].strip(' '))
                    if not level in level_names:
                        print("Unknown level! "+str(level))
                        print(name)
                        #import pdb; pdb.set_trace()
                    else:
                        self.stash_levels[scode]=level
        return()


        
    def read_stash_mappings(self):
        self.cf_to_stash = configparser.ConfigParser()   
        #turn of the lowercaseization of keys
        self.cf_to_stash.optionxform = lambda option: option
        #configFilePath = 'common_mappings.cfg'
        configFilePath = main_config['main']['mappings']
        if ',' in configFilePath:
           #this contains multiple config files, split into a list
           configFilePaths=split(',',configFilePath)
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
                    
                self.rose[cmip6_time]=cmip6_tim_dom_full
                self.rose_time_domain_mappings[freq]=this_cmip6['tim_name']
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
                      else:#
                         if this_domain in self.use_matrix[this_time]:
                            #this domain already exist here
                            if not this_use in self.use_matrix[this_time][this_domain]:
                               #only add if this use isn't already in the list
                               self.use_matrix[this_time][this_domain].append(this_use)
                         else:
                            #domain not here yet
                            self.use_matrix[this_time][this_domain]=[this_use]
                plog(this_time+" added to Use Matrix")


            else:
                #print("Time domain found in Rose")
                plog(self.rose[rose_tim_found]['tim_name']+" found in Rose")
                self.rose_time_domain_mappings[freq]=self.rose[rose_tim_found]['tim_name']
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
        for space in self.space_mappings:
            this_space_dom=self.space_mappings[space]
            #loop over all cmip6 space domains
            space_name_found=False

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
                self.rose[cmip6_space]=cmip6_space_dom_full
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

           
    
    def add_stash(self,stash_code,time_domain_cf,spatial_domain_cf):
        #checks to see if the stash_code (e.g m01s03i236) exists in the rose config object
        #and that this stash code is output using the correct time and spatial domains
        #if not, adds the stash code and any required domains

        time_domain=self.rose_time_domain_mappings[time_domain_cf]
        spatial_domain=self.rose_space_domain_mappings[spatial_domain_cf]


        #extract all stash requests from rose
        pattern = r'm(\d{2})s(\d{2})i(\d{3})'
        model,isec,item=re.match(pattern,stash_code).groups()
        this_stash_level=self.rose_get_dom_level(spatial_domain)
        sc_level=self.stash_levels[stash_code]
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
                plog("Request is for "+stash_code+" using "+spatial_domain+" but this uses level="+str(this_stash_level)+" but "+str(stash_code)+" is output on level="+str(self.stash_levels[stash_code]))
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
            print("Need to add stash code")
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
            hex_uuid = format(int(uuid.uuid4().hex[:8], 16), 'x')

            #namelist_name="[!namelist:umstash_streq("+isec+item+"_"+hex_uuid+")]"
            namelist_name="namelist:umstash_streq("+isec+item+"_"+hex_uuid+")"
            new_stash={'dom_name':spatial_domain,
                       'isec':isec1,
                       'item':item1,
                       'package':'EXTRA',
                       'tim_name':time_domain,
                       'use_name':usage
                       }
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

        
        #self.read_rose_app_conf(file+'/'+ice_conf)
        self.rose_cice=main_config['user']['job_path']+'app/nemo_cice/rose-app.conf'
        self.cice_diagnostics_file=main_config['main']['cice_diags']
        self.rose,self.rose_header=self.read_rose_app_conf(self.rose_cice)
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


    

    def addIceDiag(self,line):
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
        

        fdiag='f_'+diag
        this_section=self.rose['namelist:icefields_nml']
        if fdiag in this_section:
           #this IS a valid CICE diagnostic

           diag_freq=this_section[fdiag]
           if not this_freq in diag_freq:
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
          
           if 'f_'+diag in self.cice_diagnostics:
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
        self.missing=[] # list of diagnostics we failed to add!
        self.added=[] # list of diagnostics we succesfully to added!

        self.ocean_diag=[]
        self.ocean_diag_filename=''
        self.read_ocean_xml()
        self.file_element_id_list=self.get_file_ids()
        #self.rose={}
        
        nemo_field_def_file=main_config['main']['nemo_def']
        if not os.path.isfile(nemo_field_def_file):
            print(nemo_field_def_file+" does not exist")
            exit()
        self.nemo_vars = ET.parse(nemo_field_def_file)
        


   
    def get_file_ids(self):
        #compile list of the file element ids
        id_list=[]
        root=self.ocean_diag.getroot()
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
        ocean_diag_file=self.rose[diag_keys[0]]['source'].split('/')[-1]
        this_dir='/'.join(self.rose_conf_file.split('/')[0:-1])
        self.ocean_xml_filename=this_dir+'/file/'+ocean_diag_file
        if not os.path.exists(self.ocean_xml_filename):
            print(self.ocean_xml_filename+' does not exist?')
            import pdb; pdb.set_trace()
        #read xml file
    
        self.ocean_diag=ET.ElementTree(file=self.ocean_xml_filename)
        #root=tree.getroot()
        self.ocean_diag_filename='app'+self.ocean_xml_filename.split('app')[-1]
        return()

        

    def addOceanDiag(self,line):

        diag=line['variable']
        freq=line['time']
        dims=line['space']

        #print(+" "+freq)

        #diag = cf variable name
        #freq= output frequency mon,day etc
        print("")

        print("-------------------------")
        print("Ocean: "+diag)

        #freq_map={'mon':'1mo', 'day':'1d'}
        this_freq=self.freq_map[freq]
        this_name_suffix=self.get_name_suffix(diag)
        if this_name_suffix==None:
            print("Couldn't find the name_suffix!")
            print("Unable to add "+diag)
            return()

        root=self.ocean_diag.getroot()
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
        root=self.ocean_diag.getroot()
        fields=root.findall(".//field[@name='"+diag+"']")
        if len(fields)==0:
            #we didn't find anything!
            #try field_def.xml
            fields=self.nemo_vars.findall(".//field[@name='"+diag+"']")
            if len(fields)==0:
                #nothing here either!
                plog(bold(diag+" is unknown cf variable!"))
                return(None)

            #Found at least one field with matching name
            #extract the name from the first field attrinb
            field_name=fields[0].attrib['name']
            field_ref=fields[0].attrib['field_ref']
            #find the diagnostics with this field id
            diags=self.nemo_vars.findall(".//field[@id='"+field_ref+"']")
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
        fields=self.nemo_vars.findall(".//field[@name='"+diag+"']")

        if len(fields)==0:
            #nothing here either!
            plog(diag+" is unknown cf variable!")
            return(None)
        #Found at least one field with matching name
        #extract the name from the first field attrinb
        field_name=fields[0].attrib['name']
        field_ref=fields[0].attrib['field_ref']

        #find the diagnostics with this field id
        diags=self.nemo_vars.findall(".//field[@id='"+field_ref+"']")
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
        new_cell_measure=deepcopy(self.ocean_diag.findall(".//variable[@name='cell_measures']")[0])
        new_diag.append(new_cell_measure)
        return(new_diag)

            

    def add_field_to_file_group(self,field,file_group,name_suffix):
        #get the <file> elements in this file_group
        root=self.ocean_diag.getroot()
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
                return()

        print("Couldn't find "+this_name_suffix+"in "+file_groups[0])
        import pdb; pdb.set_trace()


    def get_file_group(self,this_freq):
        root=self.ocean_diag.getroot()
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
        self.ocean_diag.write(output_file)
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
args = parser.parse_args()

if args.config:
    conf_file=args.config
else:
    conf_file='cf_to_um.conf'

if args.stash:
    stash_type=args.stash
else:
    stash_type='um'




#read in main config file
main_config=read_config(conf_file)


#setup logging file
start_logging()
#Now rose should have all the required time and space domains defined as in the freq_mappings and space_mappings
plog("Adding to the "+bold(stash_type)+" STASH, Nemo and CICE diagnostics")
print("------------")



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
    realm=line['realm']
    #if realm in um_realms:
    if any(element in um_realms for element in realm.split(' ')):
        #all realms here use STASH
        print(".")
        um.add_cf_diagnostic(line)
    elif 'ocean' in realm:
        #add ocean diagnostic
        nemo.addOceanDiag(line)
        #print()
    elif 'seaIce' in realm:
        #add Sea Ice diagnostic
        cice.addIceDiag(line)
        #print()
    else:
        print("Unknown Realm?")
        print(line['realm'])
        import pdb; pdb.set_trace()

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
