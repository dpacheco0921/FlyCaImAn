import numpy as np
from scipy.stats import linregress

def regress(c1_, c2_, mask, step):
    c1 = c1_[mask]
    c2 = c2_[mask]
    m,yint,_,_,_ = linregress(c1,c2)
    res = m*c1+yint - c2
    pos,neg = res>0, res<0
    pc_pos = np.percentile(np.abs(res[pos]), 100.-step)
    pc_neg = np.percentile(np.abs(res[neg]), 100.-step)
    pc_dif = pc_pos-pc_neg
    return m,yint,pc_dif,pos,neg,res,pc_pos,pc_neg

def make_mask(ch1, ch2, slope_thresh=0.4, outlier_std_thresh=15, step=0.05, despeckle_kernel=5, initial_window=100, verbose=True, debug=True):
    # ch1 represents signal channel
    c1,c2 = [i.ravel() for i in [ch1,ch2]]

    # remove major outliers
    outliers = c1>c1.mean()+outlier_std_thresh*c1.std()
    c1[outliers] = 0
    c2[outliers] = 0

    # generate initial mask and regression
    i,m,win = 0,0,initial_window
    while m<slope_thresh:
        m,_,_,_,_ = linregress(c1[(c1>=i*win) & (c1<(i+1)*win)], c2[(c1>=i*win) & (c1<(i+1)*win)])
        i += 1
    mask = (c1>=i*wi)
    m,yint,dif,pos,neg,res,pc_pos,pc_neg = regress(c1,c2,mask,step=step)

    if debug:
        return mask,m

    # regression iterations
    i = 0
    while True:
        last_mask = mask.copy()
        last_par = [m,yint]
        if dif < 0 and m > slope_thresh:
            break
        exclude_idxs = np.arange(len(mask))[mask][res > pc_pos]
        mask[exclude_idxs] = False
        m,yint,dif,pos,neg,res,pc_pos,pc_neg = regress(c1,c2,mask,step=step)
        i += 1

        if verbose:
            print(i, dif, m)

    mask = last_mask
    m,yint = last_par

    maskresid = m*c1[mask]+yint - c2[mask]
    residstd = np.std(maskresid)

    resid = (m*c1+yint - c2)
    expr = resid / residstd # expression normalized to std of residuals in this brain
    expr[(resid<0) | (mask)] = 0

    mask_im = np.zeros(ch1.shape, dtype=np.float32)
    mask_im.flat[expr>0] = expr[expr>0]
    mask_im = np.array([cv2.medianBlur(i,despeckle_kernel) for i in mask_im], dtype=np.float32)
    mask_im = mask_im[:,None,...]
    
    if debug:
        fig,axs = pl.subplots(2,1)
        axs = axs.ravel()
        axs[0].scatter(c1[::100], c2[::100], marker='x', c=expr[::100], s=1./(mask.astype(int)[::100]+1)**2)
        x = np.array([c1[mask].min(), np.percentile(c1[mask], 99)])
        axs[0].plot(x, m*x+yint, color='k')
        imsh = axs[1].imshow(mask_im[:,0,...].max(axis=0), vmax=80.0)
        fig.colorbar(imsh)
        pl.savefig('/jukebox/wang/deverett/debug/{}.png'.format(name))

    return mask_im
