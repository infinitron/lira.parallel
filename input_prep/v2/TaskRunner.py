from sherpa.models.model import SimulFitModel
from SimUtils import SimUtils
from ImUtils import ImUtils
from ConfigUtils import ConfigUtils
from SpecUtils import SpecUtils
from TypeDefs import PointCel,PointPhys,FluxObj #LiraInputConfig

class TaskRunner(ConfigUtils):
    def __init__(self, config_file: str) -> None:
        self.spectral_fit_obj = None
        self.input_image_details = None
        self.psf_image_details = None
        super().__init__(config_file)

    def create_input_image(self):
        """
        Create an input image using the configuration object
        """
        self.input_image_details = ImUtils.bin_image(**dict(**self.config,nsize=self.config['inp_size']))

    def prepare_spectra(self):
        """
        Extract and fit the spectra of the core
        """
        self.spectral_fit_obj = SpecUtils.extract_and_fit_spectra(**self.config)

    def create_psf_image(self):
        """
        Simulate a PSF of the observation
        """
        psf_events_file = SimUtils.simulate_chandra_psf(**self.config,flux_obj = self.spectral_fit_obj)

        self.psf_image_detailsImUtils.bin_image(**dict(**self.config_obj,evt_file=psf_events_file),psf=True,nsize=self.config['psf_size'])
    
    def simulate_baseline_and_replicates(self):
        SimUtils.simulate_null_images(self.input_image_details,self.psf_image_details,**self.config_obj)

    def create_LIRA_inputs(self):
        self.create_input_image()
        self.prepare_spectra()
        self.create_psf_image()
        self.simulate_baseline_and_replicates()



    


