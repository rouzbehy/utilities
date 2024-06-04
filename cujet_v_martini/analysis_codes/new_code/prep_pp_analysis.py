import os

#systems = ['pp_200','pp_2760','pp_5020']
systems = ['pp_2760']

##charged_hadrons.py  jet_shape.py    observables.py 
##dijet_asymmetry.py  jet_frag_func.py  jet_spectra.py  

with open("run_list_pp.txt", 'w') as f:
    for syst in systems:
        for maxT in ['20_highstat']:#[20, 100, 150, 200]:
            #f.write(f"python3 charged_hadrons.py {syst} {maxT}\n")
            f.write(f"python3 jet_spectra.py {syst} {maxT}\n")
            #f.write(f"python3 jet_shape.py {syst} {maxT}\n")
            #f.write(f"python3 jet_frag_func.py {syst} {maxT}\n")
            #f.write(f"python3 dijet_asymmetry.py {syst} {maxT}\n")
            #f.write(f"python3 photon_spec.py {syst} {maxT}\n")

