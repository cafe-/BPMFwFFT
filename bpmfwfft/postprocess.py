"""
"""
import pickle

import numpy as np
import netCDF4

from Md_OpenMM import openmm_energy
from Md_sander import sander_energy
import IO 

#KB = 8.3144621E-3/4.184
KB = 0.001987204134799235       # kcal/mol/K

def cal_pot_energy(prmtop_file, crd, phase, tmp_dir):
    """
    combine openmm_energy() and sander_energy()
    """
    phase_prefix = phase.split("_")[0]
    if phase_prefix == "OpenMM":
        return openmm_energy(prmtop_file, crd, phase)
    elif phase_prefix == "sander":
        return sander_energy(prmtop_file, crd, phase, tmp_dir, no_bonded=True)
    else:
        raise RuntimeError("Unknown phase %s"%phase)

def bootstrapping(data, nrepetitions):
    """
    data    :   1D ndarray of float
    nrepetitions    :   int, number of repetitions 

    return
        sums    :   1D ndarray of float of shape (repetitions)
                :   sums of ramdom samples drawn from data with replacement
    """
    assert len(data.shape) == 1, "data must be 1D array"
    data = np.array(data, dtype=float)
    sample_size = data.shape[0]
    sums = []
    for i in range(nrepetitions):
        sel_ind = np.random.choice(sample_size, size=sample_size, replace=True)
        sums.append(data[sel_ind].sum())
    return np.array(sums, dtype=float)


class PostProcess:
    def  __init__(self, rec_prmtop, lig_prmtop, complex_prmtop,
                        sampling_nc_file, 
                        solvent_phases,
                        nr_resampled_complexes,
                        randomly_translate_complex,
                        temperature,
                        sander_tmp_dir):
        """
        rec_prmtop  :   str, name of receptor prmtop file
        lig_prmtop  :   str, name of ligand prmtop file
        complex_prmtop  :   str, name of complex prmtop file
        sampling_nc_file    :   str, name of nc file
        solvent_phases  :   list of str, names of OpenMM solvent phases
        nr_resampled_complexes  :   int
        randomly_translate_complex  :   bool      # TODO may be removed, not useful
        temperature :   float
        sander_tmp_dir  :   str, output dir, needed to put temp files to run sander
        """
        self._rec_prmtop = rec_prmtop
        self._lig_prmtop = lig_prmtop
        self._complex_prmtop = complex_prmtop
        self._sander_tmp_dir = sander_tmp_dir

        self._solvent_phases = solvent_phases
        # get gas phases corresponding to the solvent phases
        self._gas_phases     = self._get_gas_phases()
        print "Phases to process", self._gas_phases, self._solvent_phases

        self._temperature    = temperature
        # needed to estimate the standard deviation
        self._bootstrapping_nrepetitions = 100      # TODO why hard coded?

        self._nc_handle = netCDF4.Dataset(sampling_nc_file, "r")
        self._check_number_resampled_energy(nr_resampled_complexes)

        self._estimate_gas_bpmf()
        if self._no_sample:
            self._cal_rec_desolv_no_sample()
            self._cal_lig_desolv_no_sample()
            self._cal_complex_solv_no_sample()
        else:
            self._cal_rec_desolv()
            self._cal_lig_desolv()
            self._cal_complex_solv(nr_resampled_complexes, randomly_translate_complex)

        self._estimate_solvent_bpmf()

    def _corresponding_gas_phase(self, solvent_phase):
        return solvent_phase.split("_")[0] + "_Gas"

    def _get_gas_phases(self):
        gas_phases = [self._corresponding_gas_phase(phase) for phase in self._solvent_phases]
        return list(set(gas_phases))

    def _check_number_resampled_energy(self, nr_resampled_complexes):
        """
        To make sure the provided nr_resampled_complexes is consistent with the one in the nc file
        """
        nr_resampled_complexes_in_nc_file = self._nc_handle.variables["resampled_energies"].shape[1]

        if nr_resampled_complexes > nr_resampled_complexes_in_nc_file:
            raise RuntimeError("nr_resampled_complexes is larger than size of resampled_energies (%d) in nc file"%(
                nr_resampled_complexes_in_nc_file ) )
        return None

    def _cal_rec_desolv(self):
        """
        set self._rec_energies[phase] -> float
        and self._rec_desol_fe[phase] -> float
        """
        rec_conf = self._nc_handle.variables["rec_positions"][:]

        self._rec_energies = {}
        for p in self._gas_phases + self._solvent_phases:
            e = cal_pot_energy(self._rec_prmtop, rec_conf, p, self._sander_tmp_dir)
            self._rec_energies[p] = e[0]
        print "Receptor potential energies"
        print self._rec_energies

        self._rec_desol_fe = {}
        for solvent_phase in self._solvent_phases:
            corresp_gas_phase = self._corresponding_gas_phase(solvent_phase)
            self._rec_desol_fe[solvent_phase] = self._rec_energies[corresp_gas_phase] - self._rec_energies[solvent_phase]

        print "Receptor desolvation energies:"
        for p, e in self._rec_desol_fe.items():
            print p, e
        return None

    def _cal_rec_desolv_no_sample(self):
        """
        set self._rec_energies[phase] -> float
        and self._rec_desol_fe[phase] -> float
        """
        self._rec_energies = {}
        for p in self._gas_phases + self._solvent_phases:
            self._rec_energies[p] = 0.
        print "Receptor potential energies"
        print self._rec_energies

        self._rec_desol_fe = {}
        for solvent_phase in self._solvent_phases:
            self._rec_desol_fe[solvent_phase] = 0.

        print "Receptor desolvation energies:"
        for p, e in self._rec_desol_fe.items():
            print p, e
        return None

    def _cal_lig_desolv(self):
        """
        set self._lig_energies[phase] -> 1D np.array of len(lig_confs)
        and self._lig_desol_fe[phase] -> float
        """
        lig_confs = self._nc_handle.variables["lig_positions"][:]

        self._lig_energies = {}
        for p in self._gas_phases + self._solvent_phases:
            self._lig_energies[p] = cal_pot_energy(self._lig_prmtop, lig_confs, p, self._sander_tmp_dir)
        print "Ligand potential energies"
        print self._lig_energies

        beta = 1./ KB/ self._temperature
        self._lig_desol_fe     = {}
        self._lig_desol_fe_std = {}

        for solvent_phase in self._solvent_phases:
            corresp_gas_phase = self._corresponding_gas_phase(solvent_phase)

            delta_E = beta * (self._lig_energies[corresp_gas_phase] - self._lig_energies[solvent_phase])
            delta_E_max = delta_E.max()
            exp_energies = np.exp(delta_E - delta_E_max)

            exp_mean = exp_energies.mean()
            self._lig_desol_fe[solvent_phase] = KB * self._temperature * (np.log(exp_mean) + delta_E_max)

            exp_means = bootstrapping(exp_energies, self._bootstrapping_nrepetitions) / float(exp_energies.shape[0])
            self._lig_desol_fe_std[solvent_phase] = (KB * self._temperature * (np.log(exp_means) + delta_E_max)).std()

        print "Ligand desolvation energies:"
        for p in self._lig_desol_fe.keys():
            print p, self._lig_desol_fe[p], "+-", self._lig_desol_fe_std[p]
        return None

    def _cal_lig_desolv_no_sample(self):
        """
        set self._lig_energies[phase] -> 1D np.array of len(lig_confs)
        and self._lig_desol_fe[phase] -> float
        """
        nr_lig_confs = self._nc_handle.variables["lig_positions"].shape[0]

        self._lig_energies = {}
        for p in self._gas_phases + self._solvent_phases:
            self._lig_energies[p] = np.zeros([nr_lig_confs], dtype=float)
        print "Ligand potential energies"
        print self._lig_energies

        self._lig_desol_fe     = {}
        self._lig_desol_fe_std = {}
        for solvent_phase in self._solvent_phases:
            self._lig_desol_fe[solvent_phase] = 0.
            self._lig_desol_fe_std[solvent_phase] = 0.

        print "Ligand desolvation energies:"
        for p in self._lig_desol_fe.keys():
            print p, self._lig_desol_fe[p], "+-", self._lig_desol_fe_std[p]
        return None

    def _resample_bound_state(self, nr_resampled_complexes):
        """
        """
        beta = 1./ KB/ self._temperature

        strata_weights  = self._nc_handle.variables["exponential_sums"][:]
        log_of_divisors = self._nc_handle.variables["log_of_divisors"][:]
        common_divisor  = log_of_divisors.max()
        strata_weights *= np.exp(log_of_divisors - common_divisor)
        strata_weights /= strata_weights.sum()

        conf_trans_inds = []
        for conf_ind in range(len(strata_weights)):
            stratum_size = int(round(strata_weights[conf_ind] * nr_resampled_complexes))
            if stratum_size > 0:
                weights = self._nc_handle.variables["resampled_energies"][conf_ind]
                e_min   = weights.min()
                weights = np.exp(-beta * (weights - e_min))
                weights /= weights.sum()
                selected_trans_ind = np.random.choice(weights.shape[0], size=stratum_size, p=weights, replace=False)

                for trans_ind in selected_trans_ind:
                    conf_trans_inds.append((conf_ind, trans_ind))
        print "conf_trans_inds", conf_trans_inds
        return conf_trans_inds

    def _translate_ligand(self, conf_ind, trans_ind):
        """
        move ligand originally at lig_com to trans_corner
        """
        lig_conf     = self._nc_handle.variables["lig_positions"][conf_ind]
        trans_corner = self._nc_handle.variables["resampled_trans_vectors"][conf_ind, trans_ind]

        i, j, k = trans_corner
        x = self._nc_handle.variables["x"][i]
        y = self._nc_handle.variables["y"][j]
        z = self._nc_handle.variables["z"][k]
        
        displacement = np.array([x, y, z], dtype=float)
        
        for atom_ind in range(len(lig_conf)):
            lig_conf[atom_ind] += displacement
        return np.array(lig_conf, dtype=float)

    def _construct_complexes(self, nr_resampled_complexes, randomly_translate_complex):
        """
        set self._resampled_conf_trans_inds: list of tuples
        and self._resampled_holo_lig_confs
        receptor + ligand
        """
        self._resampled_conf_trans_inds = self._resample_bound_state(nr_resampled_complexes)
        rec_crd = self._nc_handle.variables["rec_positions"][:]

        complexes = []
        self._resampled_holo_lig_confs = []

        for conf_ind, trans_ind in self._resampled_conf_trans_inds:
            lig_crd = self._translate_ligand(conf_ind, trans_ind)
            self._resampled_holo_lig_confs.append(lig_crd)

            complexes.append(np.array(list(rec_crd) + list(lig_crd), dtype=float))

        self._resampled_holo_lig_confs = np.array(self._resampled_holo_lig_confs, dtype=float)

        if randomly_translate_complex:
            print "Randomly translate complexes"
            for i in range(len(complexes)):
                complexes[i] += np.random.rand(1)[0] / 2.

        return np.array(complexes, dtype=float)

    def _cal_complex_solv(self, nr_resampled_complexes, randomly_translate_complex):
        """
        set self._mean_interaction_energies[phase] -> float
        and self._min_interaction_energies[phase] -> float
        and self._complex_sol_fe
        """
        complex_confs = self._construct_complexes(nr_resampled_complexes, randomly_translate_complex)

        complex_energies = {}
        for p in self._gas_phases + self._solvent_phases:
            complex_energies[p] = cal_pot_energy(self._complex_prmtop, complex_confs, p, self._sander_tmp_dir)
        print "Complex potential energies"
        print complex_energies

        self._mean_interaction_energies = {}
        self._min_interaction_energies  = {}
        self._std_interaction_energies  = {}

        for p in self._gas_phases + self._solvent_phases:
            inter_energies = np.zeros(len(complex_energies[p]), dtype=float)

            for i in range(len(complex_energies[p])):
                conf_ind = self._resampled_conf_trans_inds[i][0]
                inter_energies[i] = complex_energies[p][i] - self._lig_energies[p][conf_ind] - self._rec_energies[p]

            self._mean_interaction_energies[p] = inter_energies.mean()
            self._std_interaction_energies[p]  = inter_energies.std()
            self._min_interaction_energies[p]  = inter_energies.min()

        beta = 1./ KB/ self._temperature
        self._complex_sol_fe = {}
        self._complex_sol_fe_std = {}

        for solvent_phase in self._solvent_phases:
            corresp_gas_phase = self._corresponding_gas_phase(solvent_phase)

            delta_E = -beta * (complex_energies[solvent_phase] - complex_energies[corresp_gas_phase])
            delta_E_max = delta_E.max()
            exp_energies = np.exp(delta_E - delta_E_max)
            exp_mean = exp_energies.mean()
            self._complex_sol_fe[solvent_phase] = -KB * self._temperature * (np.log(exp_mean) + delta_E_max)

            exp_means = bootstrapping(exp_energies, self._bootstrapping_nrepetitions) / float(exp_energies.shape[0])
            self._complex_sol_fe_std[solvent_phase] = (-KB * self._temperature * (np.log(exp_means) + delta_E_max)).std()
        
        print "Complex solvation energies:"
        for p in self._complex_sol_fe.keys():
            print p, self._complex_sol_fe[p], "+-", self._complex_sol_fe_std[p]
        return None

    def _cal_complex_solv_no_sample(self):
        """
        set self._mean_interaction_energies[phase] -> float
        and self._min_interaction_energies[phase] -> float
        and self._complex_sol_fe
        """
        self._resampled_conf_trans_inds = []
        self._resampled_holo_lig_confs = []

        self._mean_interaction_energies = {}
        self._min_interaction_energies  = {}
        self._std_interaction_energies  = {}

        for p in self._gas_phases + self._solvent_phases:
            self._mean_interaction_energies[p] = np.inf
            self._std_interaction_energies[p]  = 0.
            self._min_interaction_energies[p]  = np.inf

        self._complex_sol_fe = {}
        self._complex_sol_fe_std = {}
        for solvent_phase in self._solvent_phases:
            self._complex_sol_fe[solvent_phase] = 0.
            self._complex_sol_fe_std[solvent_phase] = 0.

        print "Complex solvation energies:"
        for p in self._complex_sol_fe.keys():
            print p, self._complex_sol_fe[p], "+-", self._complex_sol_fe_std[p]
        return None

    def _estimate_gas_bpmf(self):
        """
        set self._bpmf[p] for p in self._gas_phases
        set self._no_sample -> bool

        TODO: read this to learn how to get gas_bpmf
        """
        v_0 = 1661.
        beta = 1./ KB/ self._temperature
        v_binding = self._nc_handle.variables["volume"][:].mean()
        print "v_binding %f"%v_binding
        correction = -KB * self._temperature * np.log(v_binding / v_0 / 8 / np.pi**2)
        print "Volume correction %f"%correction

        nr_grid_points = np.array(self._nc_handle.variables["nr_grid_points"][:], dtype=float)
        number_of_samples = nr_grid_points.sum()

        exponential_sums = self._nc_handle.variables["exponential_sums"][:]
        log_of_divisors  = self._nc_handle.variables["log_of_divisors"][:]
        common_divisor  = log_of_divisors.max()

        exponential_sums *= np.exp(log_of_divisors - common_divisor)
        exp_mean = exponential_sums.sum() / number_of_samples

        gas_bpmf = -KB * self._temperature * (np.log(exp_mean) + common_divisor)
        gas_bpmf += correction

        # TODO check number_of_samples
        exp_means = bootstrapping(exponential_sums, self._bootstrapping_nrepetitions) / number_of_samples
        self._gas_bpmf_std = (-KB * self._temperature * (np.log(exp_means) + common_divisor)).std()
        if str(self._gas_bpmf_std) == "nan":
            self._gas_bpmf_std = 0.

        self._bpmf = {}
        for p in self._gas_phases:
            self._bpmf[p] = gas_bpmf

        self._no_sample = False
        if gas_bpmf == np.inf:
            self._no_sample = True
        return None

    def _estimate_solvent_bpmf(self):
        """
        set self._bpmf[p] for p in self._solvent_phases
        """
        for solvent_phase in self._solvent_phases:
            corresp_gas_phase = self._corresponding_gas_phase(solvent_phase)
            self._bpmf[solvent_phase] = self._bpmf[corresp_gas_phase] + self._rec_desol_fe[solvent_phase] + \
                                        self._lig_desol_fe[solvent_phase] + self._complex_sol_fe[solvent_phase]
        print "Binding potential of mean force"
        for p, e in self._bpmf.items():
            print p, e
        return None

    def get_bpmf(self):
        return self._bpmf

    def pickle_bpmf(self, file):
        data = {"bpmf" : self._bpmf, 
                "mean_Psi" : self._mean_interaction_energies,
                "min_Psi" : self._min_interaction_energies}

        standard_dev = {}
        standard_dev["interactions_energies"] = self._std_interaction_energies
        standard_dev["gas_bpmf"] = self._gas_bpmf_std
        standard_dev["lig_desolvation"] = self._lig_desol_fe_std
        standard_dev["complex_solvation"] = self._complex_sol_fe_std

        data["std"] = standard_dev

        pickle.dump(data, open(file, "w"))
        return None

    def write_rececptor_pdb(self, file):
        rec_crd = self._nc_handle.variables["rec_positions"][:]
        IO.write_pdb(self._rec_prmtop, rec_crd, file, "w")
        return

    def write_resampled_ligand_pdb(self, file):
        open(file, "w")
        for conf in self._resampled_holo_lig_confs:
            IO.write_pdb(self._lig_prmtop, conf, file, "a")
        return None


class PostProcess_PL(PostProcess):
    def _check_number_resampled_energy(self, nr_resampled_complexes):
        return None

    def _resample_bound_state(self, nr_resampled_complexes):
        """
        TODO
        """
        beta = 1./ KB/ self._temperature

        strata_weights  = self._nc_handle.variables["exponential_sums"][:]
        log_of_divisors = self._nc_handle.variables["log_of_divisors"][:]
        common_divisor  = log_of_divisors.max()
        strata_weights *= np.exp(log_of_divisors - common_divisor)
        if strata_weights.sum() != 0:
            strata_weights /= strata_weights.sum()

        conf_trans_inds = []
        for conf_ind in range(len(strata_weights)):

            key = "resampled_energies_%d"%conf_ind
            if key in self._nc_handle.variables.keys():

                stratum_size = int(round(strata_weights[conf_ind] * nr_resampled_complexes))
                if stratum_size > 0:

                    weights = self._nc_handle.variables[key][:]
                    e_min   = weights.min()
                    weights = np.exp(-beta * (weights - e_min))
                    weights /= weights.sum()

                    sample_size = min(weights.shape[0], stratum_size)
                    selected_trans_ind = np.random.choice(weights.shape[0], size=sample_size, p=weights, replace=False)
                    if sample_size < stratum_size:
                        residual_weight = float(stratum_size - sample_size)
                        residual_weight /= nr_resampled_complexes

                        if conf_ind+1 < len(strata_weights):
                            strata_weights[conf_ind+1:] += residual_weight / len(strata_weights[conf_ind+1:])
                        strata_weights[conf_ind] = float(sample_size) / nr_resampled_complexes

                    for trans_ind in selected_trans_ind:
                        conf_trans_inds.append((conf_ind, trans_ind))
                else:
                    if conf_ind+1 < len(strata_weights):
                        strata_weights[conf_ind+1:] += strata_weights[conf_ind] / len(strata_weights[conf_ind+1:])
                        strata_weights[conf_ind] = 0. 

                strata_weights /= strata_weights.sum()

        conf_trans_inds = conf_trans_inds[:nr_resampled_complexes]
        print "number of resampled complexes", len(conf_trans_inds)
        print "conf_trans_inds", conf_trans_inds
        return conf_trans_inds

    def _translate_ligand(self, conf_ind, trans_ind):
        """
        read from different key for the translational displacement
        """
        lig_conf     = self._nc_handle.variables["lig_positions"][conf_ind]

        key = "resampled_trans_vectors_%d"%conf_ind
        trans_corner = self._nc_handle.variables[key][trans_ind]

        i, j, k = trans_corner
        x = self._nc_handle.variables["x"][i]
        y = self._nc_handle.variables["y"][j]
        z = self._nc_handle.variables["z"][k]

        displacement = np.array([x, y, z], dtype=float)

        for atom_ind in range(len(lig_conf)):
            lig_conf[atom_ind] += displacement
        return np.array(lig_conf, dtype=float)

    def _cal_complex_solv(self, nr_resampled_complexes, randomly_translate_complex):
        """
        deal with case when there is no reconstructed complex
        """

        self._mean_interaction_energies = {}
        self._min_interaction_energies  = {}
        self._std_interaction_energies  = {}

        self._complex_sol_fe = {}
        self._complex_sol_fe_std = {}

        complex_confs = self._construct_complexes(nr_resampled_complexes, randomly_translate_complex)
        if complex_confs.shape[0] == 0:

            for p in self._gas_phases + self._solvent_phases:
                self._mean_interaction_energies[p] = np.inf
                self._min_interaction_energies[p]  = np.inf
                self._std_interaction_energies[p]  = 0.

            for solvent_phase in self._solvent_phases:
                self._complex_sol_fe[solvent_phase] = np.inf
                self._complex_sol_fe_std[solvent_phase] = 0.
        else:

            complex_energies = {}
            for p in self._gas_phases + self._solvent_phases:
                complex_energies[p] = cal_pot_energy(self._complex_prmtop, complex_confs, p, self._sander_tmp_dir)
            print "Complex potential energies", complex_energies

            for p in self._gas_phases + self._solvent_phases:
                inter_energies = np.zeros(len(complex_energies[p]), dtype=float)

                for i in range(len(complex_energies[p])):
                    conf_ind = self._resampled_conf_trans_inds[i][0]
                    inter_energies[i] = complex_energies[p][i] - self._lig_energies[p][conf_ind] - self._rec_energies[p]

                self._mean_interaction_energies[p] = inter_energies.mean()
                self._std_interaction_energies[p]  = inter_energies.std()
                self._min_interaction_energies[p]  = inter_energies.min()

            beta = 1./ KB/ self._temperature
            for solvent_phase in self._solvent_phases:
                corresp_gas_phase = self._corresponding_gas_phase(solvent_phase)

                delta_E = -beta * (complex_energies[solvent_phase] - complex_energies[corresp_gas_phase])
                delta_E_max = delta_E.max()
                exp_energies = np.exp(delta_E - delta_E_max)
                exp_mean = exp_energies.mean()
                self._complex_sol_fe[solvent_phase] = -KB * self._temperature * (np.log(exp_mean) + delta_E_max)

                exp_means = bootstrapping(exp_energies, self._bootstrapping_nrepetitions) / float(exp_energies.shape[0])
                self._complex_sol_fe_std[solvent_phase] = (-KB * self._temperature * (np.log(exp_means) + delta_E_max)).std()

        print "Complex solvation energies:"
        for p in self._complex_sol_fe.keys():
            print p, self._complex_sol_fe[p], "+-", self._complex_sol_fe_std[p]
        return None

#-----

