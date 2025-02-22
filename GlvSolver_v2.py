import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from scipy.interpolate import CubicSpline
from copy import deepcopy
import numpy as np
import arviz as az

###############################
# Function: compute d/dt log(y)
###############################
def compute_dlogy_dt(
    df_input,                 # dataframe: columns inlcude SubjectID, SampleID, Timepoint, Perturbation1-n, Taxon1-n
    taxa2include=None,        # list: taxa to be simulated in the gLV model
    numtaxa2include=None,     # int: number of most abundant taxa to be simulated in the gLV model (this option is only active when taxa2include is None)
    method='spline'           # str: available methods include spline and gradient (central difference)
):
    assert len(df_input)>0
    assert len(set(df_input.SampleID))==len(df_input) # each SampleID should be unique

    # only one sample per time point for each subject
    for sid in set(df_input.SubjectID):
        timepoints = list(df_input[df_input.SubjectID==sid].sort_values(by='Timepoint').Timepoint)
        assert pd.Series(timepoints).is_unique
        assert pd.Series(timepoints).is_monotonic_increasing

    # Remove time series that have a single sample (at least two time points are required)
    df_input2 = deepcopy(df_input)
    subjects_to_remove = []
    for sid in set(df_input2.SubjectID):
        if len(set(df_input2[df_input2.SubjectID==sid].SampleID)) == 1:
            subjects_to_remove.append(sid)
    if len(subjects_to_remove)>0:
        df_input2 = df_input2[~df_input2.SubjectID.isin(subjects_to_remove)]
        assert len(df_input2)>0

    # Select taxa to be simulated in the gLV model
    taxa_columns = [col for col in df_input2.columns if col.startswith('Taxon_')]
    other_columns = [col for col in df_input2.columns if not col.startswith('Taxon_')]
    if taxa2include is None:
        if numtaxa2include is not None:
            assert numtaxa2include>=1
            df_abun2 = deepcopy(df_input2[taxa_columns].T) # taxa by samples
            df_abun2['mean'] = df_abun2.mean(axis=1)
            df_abun2 = df_abun2.sort_values(by=['mean'],axis=0,ascending=False)
            taxa_included = sorted(list(df_abun2.index)[0:numtaxa2include])
            df_input2 = df_input2[other_columns+taxa_included]
        else:
            taxa_included = taxa_columns # include all taxa
    else:
        assert len(taxa2include)>0
        taxa_included = sorted(list(set(taxa2include).intersection(set(taxa_columns))))
        num_taxa_included = len(taxa_included)
        assert num_taxa_included>0
        df_input2 = df_input2[other_columns+taxa_included]

    # Glv is supposed to run on absolute abundance
    # Normalize abundance by the maximum value
    scaling_factor = df_input2[taxa_included].max().max()
    df_input2[taxa_included] = df_input2[taxa_included]/scaling_factor

    # Add pseudo abundance: replace 0 with the minimum abundance within each sample
    df_input2 = df_input2.set_index('SampleID')
    for sid in set(df_input2.index):
        min_abun = np.min([df_input2.loc[sid,taxon] for taxon in taxa_included if df_input2.loc[sid,taxon]>0])
        assert min_abun>0
        taxa_w_zero_abun = [taxon for taxon in taxa_included if df_input2.loc[sid,taxon]==0]
        df_input2.loc[sid,taxa_w_zero_abun] = min_abun
    assert np.count_nonzero(df_input2[taxa_included]) == df_input2[taxa_included].shape[0]*df_input2[taxa_included].shape[1]

    # Apply log transform
    df_y = deepcopy(df_input2[taxa_included])
    df_logy = deepcopy(np.log(df_y))
    dlogydt_columns = ['DLOGDT_'+col for col in df_y.columns]
    df_dlogydt = pd.DataFrame(index=df_y.index, columns=dlogydt_columns)

    # Calculate log-derivatives (dlogy/dt)
    for subid in list(set(df_input2.SubjectID)):
        df_tmp = df_input2[df_input2.SubjectID==subid].sort_values('Timepoint')
        timepoints = list(df_tmp.Timepoint)
        samples = list(df_tmp.index)
        for taxon in taxa_included:
            if method.lower() == 'spline':
                cs = CubicSpline(timepoints, list(df_logy.loc[samples,taxon]))
                df_dlogydt.loc[samples,'DLOGDT_'+taxon] = cs(timepoints, 1)
            elif method.lower() == 'gradient':
                df_dlogydt.loc[samples,'DLOGDT_'+taxon] = np.gradient(list(df_logy.loc[samples,taxon]), timepoints)
            else:
                raise Exception("%s is an unknown method for computing log-derivatives."%method)

    # add dlogydt to the input table
    df_output = pd.merge(df_input2, df_dlogydt, left_index=True, right_index=True, how='inner').reset_index()
    df_output = df_output[other_columns+taxa_included+dlogydt_columns]
    return df_output, scaling_factor


#######################################
# Function: Generate regression matrics
#######################################
def generate_XY_matrics(
    df_input # dataframe: columns includ SubjectID, SampleID, Timepoint, Perturbation1-n, Taxon1-n, DLOGDT_Taxon1-n
):
    # get sample ids, taxa ids, perturbation ids
    sample_ids = list(df_input.index)
    num_samples = len(sample_ids)
    taxa_ids = list([col for col in df_input.columns if col.startswith('Taxon_')])
    num_taxa = len(taxa_ids)
    perturbation_ids = [col for col in df_input.columns if col.startswith('Perturbation_')]
    num_perturbations = len(perturbation_ids)

    # get matrix of log derivative
    df_dlogydt = df_input.loc[sample_ids, ['DLOGDT_'+tid for tid in taxa_ids]]
    assert df_dlogydt.isna().sum().sum()==0 # make sure no missing values

    # get absolute abundance
    df_y = df_input.loc[sample_ids, taxa_ids]
    assert df_y.isna().sum().sum()==0 # make sure no missing values

    # get strength of perturbation
    df_perturbation = df_input.loc[sample_ids, perturbation_ids]

    # create Y matrix
    Ymat = df_dlogydt.values
    Ymat = Ymat.flatten(order='F') # column-major order (i.e., first column, then second column, then third column, and so on)

    # create X matrix
    Xmat = np.zeros(shape=(num_taxa*num_samples, (1+num_taxa+num_perturbations)*num_taxa))
    for k in np.arange(num_taxa):
        start_index = k*(1+num_taxa+num_perturbations)
        Xmat[k*num_samples:(k+1)*num_samples,start_index] = 1.0 # basal growth rate
        Xmat[k*num_samples:(k+1)*num_samples,start_index+1:start_index+1+num_taxa] = df_y.values # pairwise interactions
        Xmat[k*num_samples:(k+1)*num_samples,start_index+1+num_taxa:start_index+1+num_taxa+num_perturbations] = df_perturbation.values # perturbations

    return Xmat, Ymat, sample_ids, taxa_ids, perturbation_ids

#############################################
# Function: Generate input files for CMD stan
#############################################
def write_stan_input_file(
    prefix,                     # str: prefix of file names
    stan_path_dir,              # str: directory where stan input files are directed
    Xmat,                       # nd array: X matrix
    Ymat,                       # nd array: Y matrix
    taxa_ids,                   # list: taxa ids
    perturbation_ids,           # list: perturbations
    pairs_to_exclude=[],        # list: interactions to be excluded from the model
    sigma_ub=100,               # float: upper bound for sigma
    normal_prior_std=100,       # float: normal prior std for alpha, beta and epsilon
    neg_self_int=False          # boolean: whether to implement negative constraint on self-self interactions
):
    # number of taxa
    num_taxa = len(taxa_ids)
    num_perturbations = len(perturbation_ids)

    #---------------------------------
    # write data to stan program files
    #---------------------------------
    json_str = '{\n"N" : %d,\n' % (len(Ymat))
    json_str += '\"dlogX\" : [%s],\n' % (','.join(list(Ymat.astype(str))))

    # growth, interaction, and perturbation matrices
    json_str += '\"growth\" : [\n'
    for i in range(Xmat.shape[0]):  # Loop over each row of Xmat
        row_values = []
        for k1 in range(num_taxa):
            row_values.append(str(Xmat[i, k1 * (1 + num_taxa + num_perturbations)]))
        json_str += '[%s],\n' % (','.join(row_values))
    json_str = json_str[:-2] + '],\n'

    json_str += '\"interaction\" : [\n'
    for i in range(Xmat.shape[0]):
        row_values = []
        for k1 in range(num_taxa):
            for k2 in range(num_taxa):
                if (taxa_ids[k1], taxa_ids[k2]) in pairs_to_exclude:
                    row_values.append("0.0")
                else:
                    row_values.append(str(Xmat[i, k1 * (1 + num_taxa + num_perturbations) + 1 + k2]))
        json_str += '[%s],\n' % (','.join(row_values))
    json_str = json_str[:-2] + '],\n'

    json_str += '\"perturbation\" : [\n'
    for i in range(Xmat.shape[0]):
        row_values = []
        for k1 in range(num_taxa):
            for k2 in range(num_perturbations):
                if (taxa_ids[k1], perturbation_ids[k2]) in pairs_to_exclude:
                    row_values.append("0.0")
                else:
                    row_values.append(str(Xmat[i, k1 * (1 + num_taxa + num_perturbations) + 1 + num_taxa + k2]))
        json_str += '[%s],\n' % (','.join(row_values))
    json_str = json_str[:-2] + '],\n'

    json_str = json_str[:-2] + '}'
    with open("%s/%s.data.json" % (stan_path_dir, prefix), "w") as text_file:
        text_file.write("%s" % json_str)

    #-------------------
    # write stan program
    #-------------------

    # function block
    # function to convert a matrix to a vector in row-wise order
    model_str = """functions {
    vector row_wise_to_vector(matrix A) {
        int num_rows = rows(A);
        int num_cols = cols(A);
        vector[num_rows * num_cols] result;

        // Fill the result vector in row-wise order
        for (i in 1:num_rows) {
            for (j in 1:num_cols) {
                result[(i - 1) * num_cols + j] = A[i, j];
            }
        }
        return result;
    }
}
"""
    # data block
    model_str += 'data {\n'
    model_str += '\tint<lower=0> N;\n'
    model_str += '\tvector[N] dlogX;\n'
    model_str += '\tmatrix[N, %d] growth;\n' % num_taxa
    model_str += '\tmatrix[N, %d] interaction;\n' % (num_taxa * num_taxa)
    model_str += '\tmatrix[N, %d] perturbation;\n' % (num_taxa * num_perturbations) # orders -> [(t1,p1), (t1,p2), ..., (t2,p1), (t2,p2), ..., (t3,p1), (t3,p2), ...]
    model_str += '}\n'

    # parameter block
    model_str += 'parameters {\n'
    model_str += '\treal<lower=0, upper=%2.2f> sigma;\n' % (sigma_ub)
    model_str += '\tvector[%d] alpha;\n' % num_taxa
    model_str += '\tmatrix[%d, %d] beta;\n' % (num_taxa, num_taxa)
    model_str += '\tmatrix[%d, %d] epsilon;\n' % (num_taxa, num_perturbations)
    model_str += '}\n'

    # transformed parameter block
    if neg_self_int:
        transformed_parameters_str = """transformed parameters {
    matrix[%d, %d] constrained_beta;
    for (i in 1:%d) {
        for (j in 1:%d) {
            if (i == j) {
                constrained_beta[i, j] = -abs(beta[i, j]);
            } else {
                constrained_beta[i, j] = beta[i, j];
            }
        }
    }
}\n"""%(num_taxa, num_taxa, num_taxa, num_taxa)
    else:
        transformed_parameters_str = """transformed parameters {
    matrix[%d, %d] constrained_beta;
    for (i in 1:%d) {
        for (j in 1:%d) {
            constrained_beta[i, j] = beta[i, j];
        }
    }
}\n"""%(num_taxa, num_taxa, num_taxa, num_taxa)
    model_str += transformed_parameters_str

    # model block
    model_str += 'model {\n'
    model_str += '\tsigma ~ uniform(0,%2.2f);\n' % (sigma_ub)
    model_str += '\talpha ~ normal(0,%2.2f);\n' % (normal_prior_std)
    for i in range(num_taxa):
        for j in range(num_taxa):
            model_str += '\tbeta[%d, %d] ~ normal(0,%2.2f);\n' % (i + 1, j + 1, normal_prior_std)
    for i in range(num_taxa):
        for j in range(num_perturbations):
            model_str += '\tepsilon[%d, %d] ~ normal(0,%2.2f);\n' % (i + 1, j + 1, normal_prior_std)

    # to_vector is column-wise!!!
    model_str += '\tvector[N] mean = growth * alpha + interaction * row_wise_to_vector(beta) + perturbation * row_wise_to_vector(epsilon);\n'

    # final model expression
    model_str += '\tdlogX ~ normal(mean, sigma);\n'
    model_str += '}'

    # Write the Stan model to file
    with open("%s/%s.stan" % (stan_path_dir, prefix), "w") as text_file:
        text_file.write("%s" % model_str)

    return

##################################
# Function: Parse CMD stan output
##################################
def parse_stan_output(
    stan_output_path,       # list: path to output stan input files
    taxa_ids,               # list: taxa ids
    perturbation_ids,       # list: perturbation ids
    pairs_to_exclude=[],    # list: interactions to be excluded from the model
    bci_cutoff=0.90,        # float: Bayesian credible interval cutoff
    scaling_factor=1,       # float: scaling factor used in the function compute_dlogy_dt
    sig_only=False          # bool: show significant coefficients only
):
    n_files = len(stan_output_path)
    fit = az.from_cmdstan(stan_output_path)

    # process posterior distribution of each variable
    res = []
    for idx1, taxa1 in enumerate(taxa_ids):
        #   gather all variables related to taxa1
        all_vars = []
        all_vars.append(['alpha_%s'%(taxa1), taxa1, 'growth', idx1, np.NaN])
        for idx2, taxa2 in enumerate(taxa_ids):
            if (taxa1, taxa2) not in pairs_to_exclude:
                all_vars.append(['beta_%s_%s'%(taxa1,taxa2), taxa1, taxa2, idx1, idx2])
        for idx2, perturbation in enumerate(perturbation_ids):
            if (taxa1, perturbation) not in pairs_to_exclude:
                all_vars.append(['epsilon_%s_%s'%(taxa1,perturbation), taxa1, perturbation, idx1, idx2])

        #   generate statistics for posterior distributions
        for var in all_vars:
            data = []
            for i in np.arange(0,n_files):
                if var[0].startswith('alpha_'):
                    data.extend(list(fit.posterior['alpha'][i][:, var[3]].values))
                elif var[0].startswith('beta_'):
                    data.extend([e/scaling_factor for e in list(fit.posterior['constrained_beta'][i][:, var[3], var[4]].values)])
                elif var[0].startswith('epsilon_'):
                    data.extend(list(fit.posterior['epsilon'][i][:, var[3], var[4]].values))
            x0,x1 = az.hdi(np.array(data), hdi_prob=bci_cutoff)
            res.append([var[0], var[1], var[2], np.mean(data), np.std(data), x0, x1, x0*x1>0])
    df_parsed = pd.DataFrame(res, columns = ['Variable','Taxon1','Taxon2','Mean','STD','BCI_left','BCI_right','Sig'])

    if sig_only==True:
        return df_parsed[df_parsed.Sig==True]
    else:
        return df_parsed
