###please ensure the version should using the following versions
### python version
from platform import python_version
print('Version of python: ',python_version()) #3.8.3 
import pandas
print('Version of pandas: ',pandas.__version__) #1.2.1 "pip install pandas==1.2.1"
import numpy
print('Version of numpy: ',numpy.__version__) #1.18.5 "pip install numpy==1.18.5"
import sklearn
print('Version of sklearn: ',sklearn.__version__) #0.24.1 "pip install scikit-learn==0.24.1"
import xgboost
print('Version of xgboost: ',xgboost.__version__) #0.90 "pip install xgboost==0.90"

## import package
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier as XGBC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score as acc, confusion_matrix as cm, average_precision_score as ap
from sklearn.metrics import precision_recall_curve as prc
from sklearn.metrics import roc_curve as rc, roc_auc_score as auc
from sklearn.impute import SimpleImputer
from sklearn.impute import KNNImputer
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import seaborn as sns

def read_data(train_dir,test_dir):
    train_data = pd.read_excel(train_dir,sheet_name=4)
    train_data = train_data[features+label]
    test_id = pd.read_csv(test_dir)
    test_data = test_id[features]
    return train_data,test_data, test_id

def find_missing_feature(data):
    if label in data.columns.tolist():
        data = data.drop(label,axis=1)
    missing_features = data.columns[data.isnull().sum(axis=0)>0]
    missing_features = list(missing_features)
    return missing_features

def train_imputation(train_data):
    missing_columns = find_missing_feature(train_data)
    train_missing = train_data.loc[:,missing_columns]
    imp_mean = SimpleImputer()
    imp_constant = SimpleImputer(strategy="constant",fill_value=0.5)
    imp_mean = imp_mean.fit_transform(train_missing)
    imp_constant = imp_constant.fit_transform(train_missing)
    train_mean = train_data.copy()
    train_mean.loc[:,missing_columns]=imp_mean
    train_constant = train_data.copy()
    train_constant.loc[:,missing_columns]=imp_constant
    imputed_train_data = pd.concat([train_mean.iloc[:,0:5],train_constant.iloc[:,5:train_constant.shape[1]]],axis=1)
    return imputed_train_data

def train_KNNimputation(train_data):
    missing_columns = find_missing_feature(train_data)
    train_missing = train_data.loc[:,missing_columns]
    imp_knn = KNNImputer(n_neighbors=5)
    imp_knn = imp_knn.fit_transform(train_missing)
    imputed_train_data = train_data.copy()
    imputed_train_data.loc[:,missing_columns]=imp_knn
    return imputed_train_data

def test_imputation(test_data, mean_list):
    imputed_test_data = test_data.copy()
    for i in range(len(features)):
        if i < 5:
            imputed_test_data[features[i]].fillna(mean_list[i],inplace=True)
        else:
            imputed_test_data[features[i]].fillna(0.5,inplace=True)
    return imputed_test_data

def test_KNNimputation(test_data):
    imp_knn = KNNImputer(n_neighbors=5)
    imp_knn = imp_knn.fit_transform(test_data)
    imputed_test_data = test_data.copy()
    imputed_test_data.loc[:,:]=imp_knn
    return imputed_test_data    

def preprocess(train_data,test_data,method='knn'):
    # for mean and 0.5 imputation
    if method=='mean':
        imputed_train_data = train_imputation(train_data)
        mean_list = np.mean(imputed_train_data)
        imputed_test_data = test_imputation(test_data,mean_list)
    # for knn imputation
    if method=='knn':
        imputed_train_data = train_KNNimputation(train_data)
        imputed_test_data = test_KNNimputation(test_data)
    train_fea = imputed_train_data[features]
    train_label = imputed_train_data[label]
    test_fea = imputed_test_data[features]
    return train_fea, train_label, test_fea

def train_model(train_fea, train_label):
    parameters = {'max_depth':range(3,8),'n_estimators': np.arange(10,60,10),
    'learning_rate':np.arange(0.01, 1, 0.01),'subsample':(0.1, 1, 0.05)}
    clf = GridSearchCV(XGBC(random_state=30), parameters, cv=5, scoring='roc_auc')
    clf.fit(train_fea, train_label)
    best_model = clf.best_estimator_
    roc_auc2=auc(train_label,best_model.predict_proba(train_fea)[:,1])
    print('best parameters', clf.best_params_)
    print('acc in validation', clf.best_score_)
    print('training acc', best_model.score(train_fea, train_label))
    print(roc_auc2)
    print(best_model)
    return best_model

def test_model(best_model, test_fea, th=0.5):
    pro = best_model.predict_proba(test_fea)[:,1]
    predict=[]
    for i in pro.tolist():
        if i>th:
            predict.append(1)
        else:
            predict.append(0)
    df_predict = pd.DataFrame(predict)
    df_predict['prob'] = pro
    return df_predict

def model_main(train_dir,test_dir): 
## DATASET
    train_data,test_data,test_id = read_data(train_dir,test_dir)
    train_fea, train_label, test_fea = preprocess(train_data, test_data)
## MODELS
    best_model = train_model(train_fea, train_label)
    test_predict = test_model(best_model, test_fea)
    return train_fea, train_label, test_fea, test_predict, best_model, test_id

def process_clin(survival_dir,test_id):
    survival = pd.read_excel(survival_dir, header =0)
    survival = survival.set_index("Unnamed: 0")
    survival = survival.loc[ test_id["Unnamed: 0"].tolist(),:]
    survival_meth = survival.loc[(survival['MGMT promoter status'] == "Unmethylated") | 
                        (survival['MGMT promoter status'] == "Methylated"),:]
    return survival, survival_meth

def survivalplot(dataframe, OSorPFS,OSPFS_yaxis, censor_column, group_column, 
                  group_1_name, group_2_name, group_1_value, group_2_value, unit, plot_1or0, savefilename):
    
    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    
    plt.rc('font', family='Helvetica')
    kmf = KaplanMeierFitter()
    OScheck=dataframe.dropna(subset = [censor_column])
    T = OScheck[OSorPFS]
    C = OScheck[censor_column]
    fig, ax = plt.subplots(figsize=(4,3),dpi=200)

    group_1 = (OScheck[group_column] == group_1_value)##
    group_2 = (OScheck[group_column] == group_2_value)#

    kmf.fit(T[group_1], event_observed=C[group_1], label= group_1_name +" "+ str(sum(group_1)))
    kmf.plot(ax=ax, ci_show=False, show_censors=True, color = "#FF2500", linewidth=2) ##B8372F
    kmf.fit(T[group_2], event_observed=C[group_2], label= group_2_name +" "+ str(sum(group_2)))
    kmf.plot(ax=ax, ci_show=False, show_censors=True, color = "#0433FF", linewidth=2) ##3F5BED

    results = logrank_test(T[group_1], T[group_2], C[group_1], C[group_2], alpha=0.95 )
    plt.xlim(xmin=0)
    plt.ylim(0,1);
    plt.legend(frameon=False, fontsize = label_size, loc=1)
    
    p_value="{0:.6f}".format(results.p_value)
    print("{0:.6f}".format(results.p_value))
    plt.xlabel(unit,fontsize = label_size)
    plt.ylabel(OSPFS_yaxis,fontsize = label_size)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    
    if plot_1or0 == 1:
        plt.savefig(savefilename,bbox_inches="tight")
        print("plt saved to: ",savefilename)
        plt.show()
    elif plot_1or0 == 0:
        plt.show()
    return p_value    

def survivalcurve(dataframe, OSorPFS,OSPFS_yaxis, censor_column, group_column, 
                  group_1_name, group_2_name, group_1_value, group_2_value, unit, plot_1or0, savefilename):
    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size
    mpl.rcParams['figure.figsize'] = [4,3]
    plt.rc('font', family='Helvetica')
    kmf = KaplanMeierFitter()
    OScheck=dataframe.dropna(subset = [censor_column])
    T = OScheck[OSorPFS]
    C = OScheck[censor_column] 
    group_1 = (OScheck[group_column] == group_1_value)
    group_2 = (OScheck[group_column] == group_2_value)
    results = logrank_test(T[group_1], T[group_2], C[group_1], C[group_2], alpha=0.95 )
    p_value="{0:.6f}".format(results.p_value)
    return p_value    

def eval_survival(test_predict, test_id, survival, survival_meth, th=0.5):
    test_predict['ID']=test_id["Unnamed: 0"]
    test_predict = test_predict.rename(columns = {0:'predict'})
    TCGA = pd.merge(survival,test_predict,how="inner",left_on=survival.index, right_on="ID")
    new_group = []
    for i in range(TCGA.shape[0]):
        if TCGA.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA['predict_th'] = new_group
    group = "predict_th"
    p1=survivalcurve(TCGA,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    p2=survivalcurve(TCGA,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    TCGA_available_meth = pd.merge(survival_meth,test_predict,how="inner",left_on=survival_meth.index, right_on="ID")
    TCGA_available_meth.index = range(TCGA_available_meth.shape[0])
    new_group = []
    for i in range(TCGA_available_meth.shape[0]):
        if TCGA_available_meth.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA_available_meth['predict_th'] = new_group
    p3=survivalcurve(TCGA_available_meth,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    p4=survivalcurve(TCGA_available_meth,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    
    TCGA_unmeth = TCGA_available_meth[TCGA_available_meth['MGMT promoter status'] == "Unmethylated"]
    TCGA_unmeth.index = range(TCGA_unmeth.shape[0])
    new_group = []
    for i in range(TCGA_unmeth.shape[0]):
        if TCGA_unmeth.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA_unmeth['predict_th'] = new_group
    nr_unmeth=sum(TCGA_unmeth['predict_th'])
    p5=survivalcurve(TCGA_unmeth,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    p6=survivalcurve(TCGA_unmeth,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    
    TCGA_meth = TCGA_available_meth[TCGA_available_meth['MGMT promoter status'] == "Methylated"]
    TCGA_meth.index = range(TCGA_meth.shape[0])
    new_group = []
    for i in range(TCGA_meth.shape[0]):
        if TCGA_meth.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA_meth['predict_th'] = new_group
    nr_meth=sum(TCGA_meth['predict_th'])
    p7=survivalcurve(TCGA_meth,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    p8=survivalcurve(TCGA_meth,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    p1=float(p1)
    p2=float(p2)
    p3=float(p3)
    p4=float(p4)
    p5=float(p5)
    p6=float(p6)
    p7=float(p7)
    p8=float(p8)  
    return p1,p2,p3,p4,p5,p6,p7,p8,nr_unmeth,nr_meth

def draw_survival(test_predict, test_id, survival, survival_meth, th=0.5):
    test_predict['ID']=test_id["Unnamed: 0"]
    test_predict = test_predict.rename(columns = {0:'predict'})
    TCGA = pd.merge(survival,test_predict,how="inner",left_on=survival.index, right_on="ID")
    new_group = []
    for i in range(TCGA.shape[0]):
        if TCGA.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA['predict_th'] = new_group
    group = "predict_th"
    survivalplot(TCGA,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    survivalplot(TCGA,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    TCGA_available_meth = pd.merge(survival_meth,test_predict,how="inner",left_on=survival_meth.index, right_on="ID")
    TCGA_available_meth.index = range(TCGA_available_meth.shape[0])
    new_group = []
    for i in range(TCGA_available_meth.shape[0]):
        if TCGA_available_meth.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA_available_meth['predict_th'] = new_group
    survivalplot(TCGA_available_meth,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    survivalplot(TCGA_available_meth,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    
    TCGA_unmeth = TCGA_available_meth[TCGA_available_meth['MGMT promoter status'] == "Unmethylated"]
    TCGA_unmeth.index = range(TCGA_unmeth.shape[0])
    new_group = []
    for i in range(TCGA_unmeth.shape[0]):
        if TCGA_unmeth.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA_unmeth['predict_th'] = new_group
    survivalplot(TCGA_unmeth,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    survivalplot(TCGA_unmeth,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    
    TCGA_meth = TCGA_available_meth[TCGA_available_meth['MGMT promoter status'] == "Methylated"]
    TCGA_meth.index = range(TCGA_meth.shape[0])
    new_group = []
    for i in range(TCGA_meth.shape[0]):
        if TCGA_meth.loc[i,'prob']>th:
            new_group.append(1)
        else:
            new_group.append(0)
    TCGA_meth['predict_th'] = new_group
    survivalplot(TCGA_meth,"Disease Free (Months)","Progression Free Survival",'Disease Free Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    survivalplot(TCGA_meth,"Overall Survival (Months)","Overall Survival",'Overall Survival Status', group,"Resistant","Sensitive",1,0,"months",
                    0,"")
    return 0

def draw_ROC(best_model, train_fea, train_label):
    model2 = best_model
    X2 = train_fea
    y2 = train_label
    ns_fpr, ns_tpr, _ = rc(y2,model2.predict_proba(X2)[:,1])
    roc_auc2=auc(y2,model2.predict_proba(X2)[:,1])
    plt.figure(figsize=(6,6),dpi=120)
    plt.rc('font', family='Helvetica')
    plt.plot([0,1],[0,1],color='grey',linewidth=2)
    plt.plot(ns_fpr, ns_tpr, color='#EA322D',marker='o', label='All main features'+'(AUC = %0.2f)' % roc_auc2,linewidth=2)
    plt.xlabel('False Positive Rate',fontdict={'size':15})
    plt.ylabel('True Positive Rate',fontdict={'size':15})
    plt.legend(loc="lower right")
    plt.show()

def draw_multiple_ROC(model, fea, label, color='#EA322D', name='All main features'):
    model.fit(fea,label)
    fpr, tpr, _ = rc(label,model.predict_proba(fea)[:,1])
    roc_auc=auc(label,model.predict_proba(fea)[:,1])
    plt.plot([0,1],[0,1],color='grey',linestyle='--',linewidth=1)
    plt.plot(fpr, tpr, color=color,label=name+'(AUC = %0.2f)' % roc_auc,linewidth=2)

if __name__ == "__main__":
    cdir = os.getcwd()

    features=['EGR4_exp', 'PAPPA_exp', 'LRRC3_exp', 'ANXA3_exp','MGMT_exp','MGMT_methylated', 'subtype_Mesenchymal','subtype_Classical', 'subtype_Proneural', 'CDK4_amp', 'CDKN2A/B_del','PDGFRA_amp', 'MDM2_amp', 'EGFR_amp', 'PTEN_del', 'PIK3R1', 'PIK3CA','TP53', 'PDGFRA', 'RB1', 'EGFR', 'PIK3CG', 'ATRX', 'PTEN', 'NF1']
    label=['Response']

    train_dir = cdir+"/data/Supplementary File 1.xlsx"
    test_dir = cdir+"/data/TCGA_testing.csv"
    survival_dir = cdir+"/data/TCGA_clinical.xlsx"
    train_data,test_data,test_id = read_data(train_dir,test_dir)
    train_fea, train_label, test_fea = preprocess(train_data, test_data,method='knn')
    survival, survival_meth = process_clin(survival_dir,test_id)
    train_survial = pd.read_csv(train_dir,sheet_name=0)
    iter_lr=0.74;iter_ss=0.35;iter_th=0.6
    parameters = {'max_depth':range(3,8),'n_estimators': np.arange(10,60,10)}
    clf = GridSearchCV(XGBC(learning_rate=iter_lr,subsample=iter_ss,random_state=30), parameters, cv=5, scoring='roc_auc')
    clf.fit(train_fea, train_label)
    best_model = clf.best_estimator_
    test_predict = test_model(best_model, test_fea, th=iter_th)
    draw_ROC(best_model, train_fea, train_label)
    draw_survival(test_predict, test_id, survival, survival_meth, th=iter_th)
    bstb = best_model.get_booster()
    bstb.save_model('./TMZep_0905.bin')

  
