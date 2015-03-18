
import json, os
from collections import Counter
os.chdir('/Users/huguesfo/Devel/bioinfotools/svm')

filenames = {'TP': 'tp.json',
    'TN': 'tn.json'}

for strategy in ['Houdayer', 'AMG-diag', 'AMG-kreftgenetikk']:
    for key, filename in filenames.iteritems():
        
        data = json.load(open(filename, 'r'))
        count = Counter([v['predict'][strategy][0]['Effect'] for v in data])
        print '%s : %s' % (key, strategy)
        print 'lost sites = %d' % count.get('lost_site', 0)
        print 'no effect = %d' % (count.get('no_effect', 0) + count.get('de_novo', 0) )
        if key == 'TP':
            print 'Sensitivity (TPR) %.2f' % (count.get('lost_site', 0) / float(len(data)))
        else:
            print 'Specificity (TNR) %.2f' % ((count.get('no_effect', 0) + count.get('de_novo', 0)) / float(len(data)))
        print '--'
