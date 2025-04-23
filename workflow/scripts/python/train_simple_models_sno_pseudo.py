#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)  # ignore all future warnings
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.utils import shuffle, resample
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from imblearn.over_sampling import SMOTE
import pickle
import collections as coll

""" Train (fit) each model on the training set using the best
    hyperparameters found by hyperparameter_tuning_cv_simple_models. Pickle
    these fitted models (into .sav files) so that they can be reused after without the
    need to retrain them all over again."""

cols = ['score_c_norm', 'score_d_norm', 'score_c_prime_norm', 'score_d_prime_norm', 
        'sno_stability_norm', 'terminal_stem_stability_norm', 
        'terminal_stem_length_score_norm', 'len_norm']

# Get best hyperparameters per model
hyperparams_df = pd.read_csv(snakemake.input.best_hyperparameters, sep='\t')
hyperparams_df = hyperparams_df.drop('accuracy_cv', axis=1)
rs = snakemake.params.random_state

def df_to_params(df):
    """ Convert a one-line dataframe into a dict of params and their value. The
        column name corresponds to the key and the value corresponds to
        the value of that param (ex: 'max_depth': 2 where max_depth was the column
        name and 2 was the value in the df)."""
    cols = list(df.columns)
    params = {}
    for col in cols:
        value = df.loc[0, col]
        params[col] = value
    return params

hyperparams = df_to_params(hyperparams_df)


# Get training set
X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id')
X_train = X_train[cols]
y_train = pd.read_csv(snakemake.input.y_train, sep='\t').drop(columns=['gene_id'])

## Try to oversample the minority class (pseudogenes)
#print('Training set class distribution before SMOTE: ', coll.Counter(y_train.target))
#sm = SMOTE(random_state=rs)
#X_train, y_train = sm.fit_resample(X_train, y_train)
#print('Training set class distribution after SMOTE: ', coll.Counter(y_train.target))

# Try to undersample the majority class (expressed CD)
#print('Training set class distribution before undersampling: ', coll.Counter(y_train.target))
#data = pd.concat([X_train.reset_index(drop=True), y_train.reset_index(drop=True)], axis=1)
#majority_class = data[data.target == 1].reset_index(drop=True)  # expressed CD
#minority_class = data[data.target == 0].reset_index(drop=True)  # snoRNA pseudogene
#maj_downsampled = resample(majority_class,
#                    replace=False,  # sample without replacement
#                    n_samples=len(minority_class),  # to match minority class
#                    random_state=rs)
#downsampled_data = pd.concat([maj_downsampled, minority_class])
#downsampled_data = shuffle(shuffle(downsampled_data, random_state=rs), random_state=rs)
#y_train = downsampled_data[['target']]
#X_train = downsampled_data.drop('target', axis=1)
#print('Training set class distribution after undersampling: ', coll.Counter(y_train.target))


# Instantiate the model defined by the 'simple_models' wildcard using the best hyperparameters
# specific to each model (logreg, svc, rf, gbm, knn)
# We balance the training to overcome class imbalance for the logreg, svc and rf; gbm and knn can't do that
if snakemake.wildcards.simple_models == "logreg":
    model = LogisticRegression(C=hyperparams['C'], solver=hyperparams['solver'],
                                random_state=rs, max_iter=500, class_weight='balanced')
elif snakemake.wildcards.simple_models == "svc":
    model = svm.SVC(C=hyperparams['C'], degree=hyperparams['degree'],
                    gamma=hyperparams['gamma'], kernel=hyperparams['kernel'],
                    random_state=rs, class_weight='balanced')
elif snakemake.wildcards.simple_models == "rf":
    model = RandomForestClassifier(max_depth=hyperparams['max_depth'],
                min_samples_leaf=hyperparams['min_samples_leaf'],
                min_samples_split=hyperparams['min_samples_split'],
                n_estimators=hyperparams['n_estimators'], random_state=rs, 
                class_weight='balanced')
elif snakemake.wildcards.simple_models == "knn":
    model = KNeighborsClassifier(n_neighbors=hyperparams['n_neighbors'],
                weights=hyperparams['weights'],
                leaf_size=hyperparams['leaf_size'], p=hyperparams['p'])
else:
    model = GradientBoostingClassifier(loss=hyperparams['loss'],
                max_depth=hyperparams['max_depth'],
                min_samples_leaf=hyperparams['min_samples_leaf'],
                min_samples_split=hyperparams['min_samples_split'],
                n_estimators=hyperparams['n_estimators'], random_state=rs)

# Train model and save training accuracy to df
model.fit(X_train, y_train.values.ravel())
print(snakemake.wildcards.simple_models)
print(model.score(X_train, y_train.values.ravel()))

acc = {}
acc[snakemake.wildcards.simple_models+'_training_accuracy'] = model.score(X_train, y_train.values.ravel())
acc_df = pd.DataFrame(acc, index=[0])
acc_df.to_csv(snakemake.output.training_accuracy, sep='\t', index=False)

# Pickle the model as a .sav file ('wb' for write in binary)
pickle.dump(model, open(snakemake.output.pickled_trained_model, 'wb'))