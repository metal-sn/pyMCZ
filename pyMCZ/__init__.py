import os

p=os.path.abspath('./')
if not os.path.exists(p+'/output'):
    os.makedirs(p+'/output')

__all__=["mcz","metallicity","metscales","testcompleteness"]
