import numpy
import re

def check_nondecreasing(arr):
    print('arr.size: ', arr.size)
    for i in range(1, arr.size):
#        if arr[i] < arr[i-1] and not math.isclose(arr[i], arr[i-1]):
        if arr[i] < arr[i-1] and (not isclose(arr[i], arr[i-1])):
            return False

    return True

def to_hist(infer_cum_x):
    infer_hist_x = infer_cum_x[0]
    infer_hist_x = numpy.append(infer_hist_x, numpy.diff(infer_cum_x))
    return infer_hist_x

def load_vldb10_dataset_format(fpath):
    with open(fpath) as f:
        histogram = sorted(list(map(int, re.findall('\d+',f.read()))))

    return histogram

def expand_h_to_hg(sizehist):
    hg = []
    for (index, count) in enumerate(sizehist):
        if count > 0:
            hg += [index]*int(count)
    print('total number of groups, len(hg): ', len(hg))
    print('total number of groups, sizehist.cumsum()[-1]: ', sizehist.cumsum()[-1])

    hg_arr = numpy.array(hg)
    #print('hg_arr: ', hg_arr)
    hg_arr.sort()
    return hg_arr

def convert_esthc_to_hg(min_size_of_groups, est_hc):
    if min_size_of_groups == 1:
        #insert 0 at the index 0
        resize_hc = numpy.zeros(est_hc.size + 1)
        resize_hc[min_size_of_groups:] = est_hc
        print('min size is 1, resized hc: ', resize_hc)

    elif min_size_of_groups == 0:
        resize_hc = est_hc
        print('min size is 0, resize_hc: ', resize_hc)
    else:
        print('min size resize exception')
        raise

    est_h = to_hist(resize_hc)
    est_hg = expand_h_to_hg(est_h)
    return est_hg
