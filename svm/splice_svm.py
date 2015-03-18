# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 20:51:29 2015

@author: huguesfo
"""

import vcf, json
from sklearn import svm, preprocessing
import numpy as np

# ------------------------------------------------------------
def filter_consensus(vcf_input_filename, vcf_output_filename):
    '''
    Get variants in the Donor or Acceptor consensus
    but not the intronic 2bp nearest the junction
    '''

    splice_region = {
        ('ds', '+'):[-2, -1, 0, 3, 4, 5],
        ('ds', '-'):[-4, -3, -2, 1, 2, 3],
        ('as', '+'):[-15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, 1, 2],
        ('as', '-'):[-1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        }

    vcf_input = open(vcf_input_filename, 'r')
    vcf_reader = vcf.Reader(vcf_input)

    vcf_output = open(vcf_output_filename, 'w')
    vcf_writer = vcf.Writer(vcf_output, vcf_reader)

    for record in vcf_reader:
        splice_type = record.INFO['type'][0] # ds|as
        strand = record.INFO['strand'][0] # +|-
        location = int(record.INFO['location'][0])
        if location in splice_region[(splice_type, strand)]:
            vcf_writer.write_record(record)

    vcf_input.close()
    vcf_output.close()


# ------------------------------------------------------------
def run_scoring(vcf_filename, json_filename):
    REFSEQ="/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta"
    REFSEQGENE="/Users/huguesfo/Documents/DATA/b37/refSeq/refGene_131119.tab"
    prog = '/Users/huguesfo/Devel/genevar/genap/annotation/splice/splice_score.py'
    cmd = ['python', prog, '-i', vcf_filename, '-o', json_filename, '--all',
           '--refseqgene', REFSEQGENE, '--refseq', REFSEQ]
    #subprocess.call(cmd)
    print 'Run this:'
    print ' '.join(cmd)

# ------------------------------------------------------------
def make_feature_matrix(tp_json_filename, tn_json_filename, algorithms = ['MaxEntScan', 'SSFL', 'GeneSplicer', 'NNSplice', 'HSF'], use_ratio = False ):
    '''
    returns
     - an array X of size [n_samples, n_features] holding the training samples,
     - an array y of class labels (strings or integers), size [n_samples]
       0 for negative (non-pathogenic)
       1 for positive (pathogenic)
    '''


    feature_matrix = []
    label_vector = []
    for idx, filename in enumerate([tn_json_filename, tp_json_filename]):
        records = json.load(open(filename, 'r'))
        for record in records:
            pos = record['authentic']['pos']
            for wild_score in record['wild']:
                if wild_score['pos'] == pos:
                    break
            for mut_score in record['mut']:
                if mut_score['pos'] == pos:
                    break
            feature_vector = []
            for algorithm in algorithms:
                if use_ratio:
                    ws = wild_score['scores'].get(algorithm, 0.0)
                    ms = mut_score['scores'].get(algorithm, 0.0)
                    if ms!=0 and ws!=0:
                        ratio = ms/float(ws) - 1.0
                    elif ms==0 and ws!=0:
                        ratio = -1.0
                    elif ms!=0 and ws==0:
                        ratio = 1.0
                    elif ms==0 and ws==0:
                        ratio = 1.0
                    feature_vector.extend([ratio])
                else:
                    feature_vector.extend([
                        wild_score['scores'].get(algorithm, 0),
                        mut_score['scores'].get(algorithm, 0),
                        ])


            feature_matrix.append(feature_vector)
            label_vector.append(idx)

    return np.array(feature_matrix), np.array(label_vector)

# ============================================================
def stuff():
    '''
    http://kukuruku.co/hub/python/introduction-to-machine-learning-with-python-andscikit-learn

    '''

    #os.chdir('/Users/huguesfo/Devel/bioinfotools/svm')
    #vcf_input_filename = '/Users/huguesfo/Documents/DATA/Splice/TP/hgmd_splice_01.vcf'
    #vcf_output_filename = 'tp_consensus.vcf'
    #filter_consensus(vcf_input_filename, vcf_output_filename)
    tp_json_filename = 'tp_consensus.json'
    #run_scoring(vcf_output_filename, tp_json_filename)

    tn_json_filename = '../validation/splice_cited_MAF_scored.json'

    X, y = make_feature_matrix(tp_json_filename, tn_json_filename)
    # normalize the data attributes
    normalized_X = preprocessing.normalize(X)
    # standardize the data attributes
    standardized_X = preprocessing.scale(X)

    # don't overfit!!
    from sklearn.cross_validation import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(standardized_X, y, test_size=0.2)

    # compute the informativeness of features.
    from sklearn import metrics
    from sklearn.ensemble import ExtraTreesClassifier
    model = ExtraTreesClassifier()
    model.fit(X_train, y)
    # display the relative importance of each attribute
    print(model.feature_importances_)
    # [ 0.12559677  0.17176821  0.09335772  0.11713894  0.03273589  0.04534374
    # 0.08758765  0.15266255  0.08398096  0.08982757]

    # search of subsets of features in order to find the best subset
    # Recursive Feature Elimination Algorithm
    from sklearn.feature_selection import RFE
    from sklearn.linear_model import LogisticRegression
    model = LogisticRegression()
    # create the RFE model and select 4 attributes
    rfe = RFE(model, 10)
    rfe = rfe.fit(X, y)
    # summarize the selection of the attributes
    print(rfe.support_)
    print(rfe.ranking_)
    '''
    4:
    [ True  True False False False False  True  True False False]
    [1 1 5 6 2 3 1 1 7 4]
    '''

    ## LogisticRegression
    from sklearn import metrics
    from sklearn.linear_model import LogisticRegression
    model = LogisticRegression()
    model.fit(X_train, y_train)
    print(model)
    # make predictions
    expected = y_test
    predicted = model.predict(X_test)
    from sklearn.metrics import r2_score
    r2_score(expected, predicted)
    # summarize the fit of the model
    print(metrics.classification_report(expected, predicted))
    print(metrics.confusion_matrix(expected, predicted))
#    model.predict(X[15])
#    model.predict_proba(X[15])


    clf = svm.SVC(kernel='linear', C=1.0)
    clf.fit(feature_matrix, label_vector)
    #clf.predict([[10, 6, 85, 72]])

    # get support vectors
    print clf.support_vectors_
    # get indices of support vectors
    print clf.support_
    # get number of support vectors for each class
    print clf.n_support_

    # dec = clf.decision_function([[10, 6, 85, 72]])
    # clf.get_params()

    # get the separating hyperplane
    w = clf.coef_[0]
    print w

    neg_predicted = clf.predict(feature_matrix[:374])
    pos_predicted = clf.predict(feature_matrix[374:])

    # "oracle" sensitivity: 86.04% (houdayer 78.03%, AMG-diag 89.66%)
    TPR = sum(pos_predicted==1) / float(len(pos_predicted))
    print TPR
    # "oracle" specificity: 95.18% (houdayer 88.40%, AMG-diag 60.00%)
    TNR = sum(neg_predicted==0) / float(len(neg_predicted))
    print TNR

# ------------------------------------------------------------
class stats(dict):
    def __init__(self, conf_matrix):
        super(stats, self).__init__()
        self.conf_matrix = conf_matrix
        self.compute()

    def compute(self):
        "compute some usefull stats from the confusion matrix"
        TP=float(self.conf_matrix[0][0])
        FP=float(self.conf_matrix[0][1])
        FN=float(self.conf_matrix[1][0])
        TN=float(self.conf_matrix[1][1])

        self['Sensitivity'] = TP / (TP+FN) # probability of positive test given that data is positive (pathogenic)

        self['Specificity'] = TN / (FP+TN) # probability of negative test given that data is negative (non pathogenic)

        self['Accuracy'] = (TP + TN) / (TP+TN+FP+FN)

        #Positive predictive value (PPV, Precision)
        self['PPV'] = TP / (TP+FP)

        #Negative predictive value (NPV)
        self['NPV'] = TN / (FN+TN)

        #Positive likelihood ratio
        self['LR+'] = (TP / (TP+FN)) / (FP/(FP+TN))

        #Negative likelihood ratio
        self['LR-'] = (FN/(TP+FN))/(TN/(FP+TN))


# ------------------------------------------------------------
def plot_Xy(X, y, lr=None):
    import matplotlib.pyplot as plt
    # 0 is non pathogenic, 1 is pathogenic

    fig, ax = plt.subplots()
    if lr:
        line_x = np.array([0, 16])
        coef = lr.coef_[0] #[ 1.92149774, -2.81144662]
        intercept = lr.intercept_ #0.75475917
        line_y = (-intercept -coef[0] * line_x) / coef[1]

        ll = ax.plot(line_x, line_y, color='k', linestyle='-', linewidth=2)

    neg = ax.scatter([x_[0] for (x_,y_) in zip(X,y) if y_==0],
                [x_[1] for (x_,y_) in zip(X,y) if y_==0],
                color='blue', marker='.')
    pos = ax.scatter([x_[0] for (x_,y_) in zip(X,y) if y_==1],
                [x_[1] for (x_,y_) in zip(X,y) if y_==1],
                color='red', marker='x')

    ax.legend((neg, pos), ('non-pathogenic','pathogenic'))
    ax.set_title('Scores distribution')
    ax.set_xlabel('MaxEntScan-wild')
    ax.set_ylabel('MaxEntScan-mut')
    ax.set_ylim([0, 16])
    ax.set_xlim([0, 16])
    ax.grid()
    ax.set_frame_on(False)
    return plt
    #plt.show()


# ============================================================
if __name__ == '__main__':
    tp_json_filename = 'tp_consensus.json'
    tn_json_filename = '../validation/splice_cited_MAF_scored.json'
    X, y = make_feature_matrix(tp_json_filename, tn_json_filename, algorithms = ['MaxEntScan'] )

    from sklearn.linear_model import LogisticRegression
    from sklearn import svm
    from sklearn.cross_validation import KFold
    from sklearn.preprocessing import StandardScaler
    from sklearn import metrics


    # Always use shuffle=True to produce random folds. If you need non-random splits, you're probably doing something wrong.
    kf = KFold(len(y), n_folds=5, shuffle=True)

    y_pred = np.zeros(len(y)) #, dtype=y.dtype) # where we'll accumulate predictions
    lr = LogisticRegression(penalty='l2',
                            dual=False,
                            tol=0.0001,
                            C=1.0,
                            fit_intercept=True, # set to False for intercept=0
                            intercept_scaling=1,
                            class_weight=None,
                            random_state=None)
    #clf = svm.SVC(kernel='linear')

    # train_index and test_index never have the same values, test_indexes never overlap
    for train_index, test_index in kf:

        # for each iteration of the for loop we'll do a test train split
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        t = StandardScaler()
        X_train = t.fit_transform(X_train)
        lr.fit(X_train, y_train) # Train on the training data
        #clf.fit(X_train, y_train)

        X_test = t.transform(X_test)
        y_pred[test_index] = lr.predict(X_test) # Predict using the test and store in y_pred
        #y_pred[test_index] = clf.predict(X_test) # Predict using the test and store in y_pred

    print metrics.r2_score(y, y_pred)
    print(metrics.classification_report(y, y_pred))
    print(metrics.confusion_matrix(y, y_pred))

    plot_Xy(X, y, lr)

    '''
    Using all 5 algos:

         precision    recall  f1-score   support

              0       0.84      0.95      0.89       374
              1       0.96      0.86      0.91       473

    avg / total       0.91      0.90      0.90       847

    [[357  17]
     [ 68 405]]

    Only MES:
             precision    recall  f1-score   support

          0       0.84      0.95      0.89       374
          1       0.95      0.86      0.90       473

avg / total       0.90      0.90      0.90       847

[[354  20]      [[TP FP]
 [ 68 405]]      [FN TN]]


    sensitivity = tp / (tp+fp) = 354 / (354.0+20) = 95%
    specificity = tn / (tn+fn) = 405 / (68+405.0) = 86%

    SVM linear: no improvement
    SVM RBF kernel: no improvement

    predict_proba(X)

    lr.coef_ : [ 1.92149774, -2.81144662]
    lr.intercept_ : 0.75475917


    clf.coef_ : [ 1.51338591, -2.01630377]
    clf.intercept_ : 0.43154236

    if interect==0:
    [[361  13]
    [ 87 386]]
    y = 0.62*x -0.00 (aka MaxEntScan -38%)
    sensitivity = tp / (tp+fp) = 361 / (361.0+13) = 96%
    specificity = tn / (tn+fn) = 386 / (87+386.0) = 81%


    '''


