import subprocess
import os,sys
if __name__ == '__main__':
    query_graphs = ['star2.g','star.g','e1.g','e2.g','e3.g','e4.g',
                    'r1.g','r2.g','r3.g','r4.g','r5.g','r6.g','r7.g','r8.g','5nodes/1.g','5nodes/12.g',
                    '5nodes/23.g','6nodes/24.g','6nodes/12.g','6nodes/1.g',
                    '7nodes/2.g','7nodes/3.g','7nodes/15.g']
    data_graphs = ['Enron.g','gowalla.g','roadNetCa.g','roadNetPa.g','roadNetTx.g']
    for d in data_graphs:
        for q in query_graphs:
            subprocess.run(['./gm','../../sm_dataset/data/ours_format/{}'.format(d),'../../sm_dataset/query/ours_format/{}'.format(q)])
