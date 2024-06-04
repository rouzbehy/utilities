template = """#!/usr/bin/env bash
#SBATCH --job-name={job_name}
#SBATCH --time=2-14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rouzbeh.modarresi-yazdi@mail.mcgill.ca
#SBATCH --output=./outs/%x_%a.out
#SBATCH --account=def-gale

cd /home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/analysis/
source venv/bin/activate
parallel --resume --joblog ./logs/{log_name} < AA_Job_List."""

systems = ['PbPb2760']
energy_loss_modules = ['cujet', 'martini']
cents = ['00-05','05-10','10-20','20-30','30-40','40-50']

for collision_syst in systems:
    for eloss in energy_loss_modules:
        for cent in cents:
            job_fname = f"{collision_syst}_{eloss}_{cent}.txt"
            log_fname = f"{collision_syst}_{eloss}_{cent}.log"
            job_name = f"{eloss}_{cent}_{collision_syst}"
            with open(f"submit_{collision_syst}_{eloss}_{cent}.sh",'w') as f:
                f.write(template.format(job_name=job_name, 
                                        log_name=log_fname,
                                        task_list_name=job_fname))
