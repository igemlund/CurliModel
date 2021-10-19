## CurliModel
This is the iGEM Lund Inhibitor model. The aim of this model was be able to evaluate different inhibitors. Using the model, one can evaluate how curli formation occurs under different bacteria concentrations and on various different time-scales. The simultations are accurate up to 10 hours of growth after which the bacteria usually start forming colonies and become harder to model. It is also possible to simulate how different inhibitors affect the curli formation. Using the modular **Inhibitor** class,  one can fit different rate parameters or investigate possible interaction mechanisms. Lastly, one can simulate several different bacteria that express different inhibiros at different distances from each other to see what production rates would be required to inhibit curli formation.

# fibrilformation.py
Includes the **FibrilFormation** classes. These are used for all simulations. There is one model that simulates uniform protein concentrations and one that includes diffusion of protein. The uniform model works great for curli at bacteria concentrations form 1000 bacteria/dm^3 and up and for smaller inhibitors. When simulating inhibitors larger than 100 kDa, or very small bacteria concentrations, one needs to include diffusion. The diffusion model is significally slower, less developped and flexible.

# Making an inhibitor
Every inhibitor needs:
- a funtion **timeStep** that describes how the inhibitor and curli concentration will be affected by the inhibitor every timestep. 
- a function **rateFunc** that describes how the inhibitor affects the curli formation rate parameters (could just return them unaffected).

# Jupyter notebooks
There are several notebooks in the bin folder that describe how we designed and used the model. These could easily be altered to work for a different protein or inhibior.





# Attributions
The work is largely based on Dr. Georg Meisl's work on the mathematics of amyloid formation. He also generously agreed to aid us implementing the CsgC chaperone. 
Reference to other work is found in the notebooks.
