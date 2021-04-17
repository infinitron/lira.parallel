from typing import Dict
from TypeDefs import LiraInputConfig
import yaml
import shutil


class ConfigUtils():
    def __init__(self,config_file:str) -> None:
        self.config=self.process_config(config_file)
    
    def process_config(self,config_file:str='lira_input_prep.yaml')->LiraInputConfig:
        """
        Parses the config file and returns the config object. The values that aren't set in the config are set to their defaults.
        :param config_file: Path to the YAML config file
        ...
        :rtype: LiraInputConfig
        """
        with open("lira_input_prep.yaml", "r") as stream:
            try:
                params = yaml.safe_load(stream)
                param_def = self.get_param_definition()
                config_object:Dict=dict()
                # go through all the params and edit the if a user supplies it
                for key, defn in param_def.items():
                    if (defn["required"]) and not key in params:
                        raise Exception("{0} is required".format(key))
                    if defn["required"] and defn["type"] != type(params[key]).__name__:
                        raise Exception("Incorrect value for {0}".format(key))
                    if key in params:
                        param_def[key]["value"] = params[key]

                #for backwards compatibility
                param_def['reg_file'] = param_def['core_reg']
                param_def['bkg_file']  = param_def['bkg_reg']
                
                for k,v in param_def.items():
                    config_object[k]=v['value']

                return LiraInputConfig(config_object)

            except yaml.YAMLError as exc:
                print(exc)

    def get_param_definition(self)->Dict:
        param_list = [
            "evt_file",
            "binsize",
            "core_reg",
            "bkg_reg",
            "n_psf_sims",
            "n_null_sims",
            "inp_size",
            "psf_size",
            "nH",
            "add_gal",
            "redshift",
            "group",
            "blur",
            "center",
            "no_core",
            "sim_baselines",
            "extract_spectra",
            "fit_spectra",
        ]
        param_types = [
            "str",
            "float",
            "str",
            "str",
            "int",
            "int",
            "int",
            "int",
            "float",
            "int",
            "float",
            "int",
            "float",
            "list",
            "bool",
            "bool",
            "bool",
            "bool",
        ]
        param_req = [
            True,
            False,
            True,
            False,
            True,
            False,
            False,
            True,
            False,
            False,
            False,
            False,
            True,
            False,
            True,
            False,
            False,
            False,
        ]
        default_params = [
            None,
            0.5,
            None,
            None,
            50,
            50,
            64,
            32,
            None,
            0.0,
            0.0,
            10,
            0.25,
            None,
            False,
            True,
            True,
            True,
        ]
        param_definition = {}
        for i in range(0, len(param_list)):
            param_definition[param_list[i]] = {
                "type": param_types[i],
                "required": param_req[i],
                "value": default_params[i],
            }

        return param_definition