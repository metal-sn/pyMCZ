import os

p=os.path.abspath('..')
if not os.path.exists(p+'\\sn_data'):
    os.makedirs(p+'\\sn_data')
mcz='MCZ'
__all__=[mcz]
