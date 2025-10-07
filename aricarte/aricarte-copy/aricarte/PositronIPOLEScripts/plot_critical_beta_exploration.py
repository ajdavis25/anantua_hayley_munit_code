import matplotlib.pyplot as plt
import numpy as np
import os
import h5py

def plotCriticalBetaTest(testFolder, output=None, xlim=[-80,80], ylim=[-80,80], magneticFieldRequirement='M'):

    if testFolder[-1] != '/':
        testFolder += '/'

    #Assuming a particular naming scheme as well as there being 9 images.
    imageFiles = os.listdir(testFolder)
    imageFiles = [image for image in imageFiles if (image.split('_')[1][0] == magneticFieldRequirement)]
    f_list = np.array([float(file.split('_')[-4]) for file in imageFiles])
    bc_list = np.array([float(file.split('_')[-3]) for file in imageFiles])
    unique_f = np.unique(f_list)
    unique_bc = np.unique(bc_list)

    fig, axarr = plt.subplots(3, 3, figsize=(8,8))
    imageMatrix = [[] for dummy in range(axarr.shape[0])]
    for f_index in range(len(unique_f)):
        for bc_index in range(len(unique_bc)):
            imageIndex = np.where((f_list == unique_f[f_index]) & (bc_list == unique_bc[bc_index]))[0][0]
            imageMatrix[f_index].append(imageFiles[int(imageIndex)])

    for row in range(axarr.shape[0]):
        for col in range(axarr.shape[1]):
            ax = axarr[row,col]

            #Open the file
            imageFile = imageMatrix[row][col]
            with h5py.File(testFolder + imageFile, 'r') as myfile:
                image = myfile['unpol'][()].transpose()
                scale = myfile['header']['scale'][()]
                dx = myfile['header']['camera']['dx'][()]
                dy = myfile['header']['camera']['dy'][()]
                dsource = myfile['header']['dsource'][()]
                Lunit = myfile['header']['units']['L_unit'][()]

            fovx_uas = dx * Lunit / dsource * 2.06265e11
            fovy_uas = dy * Lunit / dsource * 2.06265e11
            extent = [ -fovx_uas/2., fovx_uas/2., -fovy_uas/2., fovy_uas/2. ]

            image *= scale
            totalIntensity = np.sum(image)

            ax.imshow(image, origin='lower', cmap='hot', vmin=0, vmax=np.max(image), extent=extent)
            ax.text(0.07, 0.1, "I={0:3.3f} Jy".format(totalIntensity), color='white', transform=ax.transAxes)
            ax.text(0.07, 0.9, r"f={0:3.2f}, $\beta_c$={1:3.2f}".format(unique_f[row], unique_bc[col]), color='white', transform=ax.transAxes)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xticks(np.linspace(-70,70,5))
            ax.set_yticks(np.linspace(-70,70,5))
            if row < axarr.shape[0]-1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('$x \ [\mu\mathrm{as}]$', fontsize=12)
            if col > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('$y \ [\mu\mathrm{as}]$', fontsize=12)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    if output is None:
        fig.show()
    else:
        fig.savefig(output)
        plt.close(fig)

if __name__ == '__main__':
    testFolder = '../ipole_output/optimizeMunit_inc10/'
    testFolders = ['../ipole_output/optimizeMunit_inc10/', '../ipole_output/optimizeMunit_inc50/', '../ipole_output/optimizeMunit_inc90/']
    for testFolder in testFolders:
        for b in ['M', 'S']:
            output = '../plots/critical_beta_exploration/critical_beta_' + b + '_' + testFolder.split('_')[-1][:-1] + '.pdf'
            plotCriticalBetaTest(testFolder, output, magneticFieldRequirement=b)
