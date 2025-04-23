#!/usr/bin/python3
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV, RandomizedSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.utils import resample, shuffle
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import precision_score, recall_score, f1_score
from imblearn.over_sampling import SMOTE
import collections as coll

""" Tune the hyperparameters of each models (Logistic regression, Support
    vector classifier, Random Forest and Gradient boosting classifier)
    before even training them, using GridSearchCV with stratified k-fold on the
    cross-validation (X_cv) set."""
cols = ['score_c_norm', 'score_d_norm', 'score_c_prime_norm', 'score_d_prime_norm', 
        'sno_stability_norm', 'terminal_stem_stability_norm', 
        'terminal_stem_length_score_norm', 'len_norm']
X_cv = pd.read_csv(snakemake.input.X_cv, sep='\t', index_col='gene_id')
X_cv = X_cv[cols]
y_cv = pd.read_csv(snakemake.input.y_cv, sep='\t').drop(columns=['gene_id'])
rs = snakemake.params.random_state
output = snakemake.output.best_hyperparameters


## Try to oversample the minority class (pseudogenes)
#print('Tuning set class distribution before SMOTE: ', coll.Counter(y_cv.target))
#sm = SMOTE(random_state=rs)
#X_cv, y_cv = sm.fit_resample(X_cv, y_cv)
#print('Tuning set class distribution after SMOTE: ', coll.Counter(y_cv.target))


## Try to undersample the majority class (expressed CD)
#print('Tuning set class distribution before undersampling: ', coll.Counter(y_cv.target))
#data = pd.concat([X_cv.reset_index(drop=True), y_cv.reset_index(drop=True)], axis=1)
#print(data)
#majority_class = data[data.target == 1].reset_index(drop=True)  # expressed CD
#minority_class = data[data.target == 0].reset_index(drop=True)  # snoRNA pseudogene
#maj_downsampled = resample(majority_class,
#                    replace=False,  # sample without replacement
#                    n_samples=len(minority_class),  # to match minority class
#                    random_state=rs)
#downsampled_data = pd.concat([maj_downsampled, minority_class])
#downsampled_data = shuffle(shuffle(downsampled_data, random_state=rs), random_state=rs)
#y_cv = downsampled_data[['target']]
#X_cv = downsampled_data.drop('target', axis=1)
#print(X_cv)
#print('Tuning set class distribution after undersampling: ', coll.Counter(y_cv.target))


# Instantiate the model defined by the 'models' wildcard
if snakemake.wildcards.simple_models == "logreg":
    model = LogisticRegression(random_state=rs, max_iter=500)
elif snakemake.wildcards.simple_models == "svc":
    model = svm.SVC(random_state=rs)
elif snakemake.wildcards.simple_models == "rf":
    model = RandomForestClassifier(random_state=rs)
elif snakemake.wildcards.simple_models == "knn":
    model = KNeighborsClassifier()
else:
    model = GradientBoostingClassifier(random_state=rs)

# Configure the cross-validation strategy (StratifiedKFold where k=3)
cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=rs)

# Define search space (i.e. in which range of params the GridSearch happens) per model
space = snakemake.params.hyperparameter_space

# Execute the gridsearch per model on the CV set
search = GridSearchCV(estimator=model, param_grid=space,
                        cv=cv, scoring="f1_macro")
search.fit(X_cv, y_cv.values.ravel())
print(snakemake.wildcards.simple_models)
print(search.best_score_)
print(search.best_params_)

# Return the best hyperparameters found by GridSearchCV, and the accuracy of each model
# fitted on the CV set with these hyperparameters into a dataframe
params_dict = search.best_params_
params_dict['accuracy_cv'] = search.best_score_
params_df = pd.DataFrame(params_dict, index=[0])
params_df.to_csv(output, sep='\t', index=False)
