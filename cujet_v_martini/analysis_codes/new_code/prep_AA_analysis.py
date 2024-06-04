import os

collision_syst = 'PbPb2760'
energy_loss_modules = ['cujet', 'martini']
cents = ['00-05','05-10','10-20','20-30','30-40','40-50']

with open("AA_Job_LIST.txt","w") as f:
    for eloss in energy_loss_modules:
        for cent in cents:
           inputs = f"{collision_syst} {eloss} {cent}"
           f.write(f"python3 charged_hadrons.py {inputs}\n")
           f.write(f"python3 jet_spectra.py {inputs}\n")
           f.write(f"python3 jet_shape.py {inputs}\n")
           #f.write(f"python3 jet_frag_func.py {inputs}\n")
           f.write(f"python3 dijet_asymmetry.py {inputs}\n")
           f.write(f"python3 photon_spec.py {inputs}\n") 
