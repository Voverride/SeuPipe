import spateo as st
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from threading import Lock
import math
import scipy.sparse
import os
from scipy.sparse import lil_matrix


def preprocess(adatasub, bin_file, tmp_dir, layer_name, prealigned, align, startx, starty, patchsize, bin_size, n_neighbor):
    #read data
    # if prealigned:
    #     adatasub = st.io.read_bgi_agg(bin_file, image_file, prealigned=True)
    # else:
    #     adatasub = st.io.read_bgi_agg(bin_file, image_file)
    if int(patchsize) > 0:
        adatasub = adatasub[int(startx):int(startx)+int(patchsize),int(starty):int(starty)+int(patchsize)].copy()
    
    layer_data = adatasub.layers[layer_name]
    if scipy.sparse.issparse(layer_data):
        adatasub.layers[layer_name] = layer_data.toarray()

    startx = str(startx)
    starty = str(starty)
    patchsize = str(patchsize)

    #print(adatasub.shape)
    # adatasub.layers['unspliced'] = adatasub.X
    patchsizex, patchsizey = adatasub.shape

    #align staining image with bins
    # before = adatasub.layers['stain'].copy()
    # if align:
    #     st.cs.refine_alignment(adatasub, mode=align, transform_layers=['stain'])

    # fig, axes = plt.subplots(ncols=2, figsize=(16, 8), tight_layout=True)
    # axes[0].imshow(before)
    # st.pl.imshow(adatasub, 'unspliced', ax=axes[0], alpha=0.6, cmap='Reds', vmax=10, use_scale=False, save_show_or_return='return')
    # axes[0].set_title('before alignment')
    # st.pl.imshow(adatasub, 'stain', ax=axes[1], use_scale=False, save_show_or_return='return')
    # st.pl.imshow(adatasub, 'unspliced', ax=axes[1], alpha=0.6, cmap='Reds', vmax=10, use_scale=False, save_show_or_return='return')
    # axes[1].set_title('after alignment')

    # try:
    #     os.mkdir('fig/')
    # except FileExistsError:
    #     print('fig folder exists')

    # plt.savefig('fig/alignment' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.png')

    #nucleus segmentation from staining image
    # fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
    # st.cs.mask_nuclei_from_stain(adatasub, otsu_classes = 4, otsu_index=1)
    # st.pl.imshow(adatasub, 'stain_mask', ax = ax)

    # plt.savefig('fig/stain_mask' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.png')

    # st.cs.find_peaks_from_mask(adatasub, 'stain', 7)
    # st.cs.watershed(adatasub, 'stain', 5, out_layer='watershed_labels')

    # fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
    # st.pl.imshow(adatasub, 'stain', save_show_or_return='return', ax=ax)
    # st.pl.imshow(adatasub, 'watershed_labels', labels=True, alpha=0.5, ax=ax)
    # plt.savefig('fig/watershed_labels' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.png')
    #adatasub.write('data/Mouse_brain_Adult_5800:8000:900:900.h5ad')
    #print(adatasub)

    adatasub.write(tmp_dir+'/spots' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.h5ad')

    print('Prepare data for neural network...')

    # watershed2x = {} # 获取每个标签对应的x坐标， {label:[x1,x2,x3...]}
    # watershed2y = {} # 获取每个标签对应的y坐标， {label:[y1,y2,y3...]}

    # for i in range(adatasub.layers['watershed_labels'].shape[0]):
    #     for j in range(adatasub.layers['watershed_labels'].shape[1]):
    #         if adatasub.layers['watershed_labels'][i, j] == 0:
    #             continue
    #         if adatasub.layers['watershed_labels'][i, j] in watershed2x:
    #             watershed2x[adatasub.layers['watershed_labels'][i, j]].append(i)
    #             watershed2y[adatasub.layers['watershed_labels'][i, j]].append(j)
    #         else:
    #             watershed2x[adatasub.layers['watershed_labels'][i, j]] = [i]
    #             watershed2y[adatasub.layers['watershed_labels'][i, j]] = [j]

    # watershed2center = {} # 获取每个标签对应的中心点， {label:[x,y]}
    # sizes = [] # 每一个label对应的像素数量列表， [size1, size2, size3...]
    # for nucleus in watershed2x:
    #     watershed2center[nucleus] = [np.mean(watershed2x[nucleus]), np.mean(watershed2y[nucleus])]
    #     sizes.append(len(watershed2x[nucleus]))

    mask_dict = defaultdict(list)
    lock = Lock()
    labels = adatasub.layers[layer_name]

    def process_cell(i, j):
        label = labels[i, j]
        if label != 0: 
            with lock:
                mask_dict[label].append((i, j))
    
    with ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(process_cell, i, j): (i, j)
            for i in range(labels.shape[0])
            for j in range(labels.shape[1])
        }

    mask2center = {}
    sizes = []
    for nucleus in mask_dict:
        array = np.array(mask_dict[nucleus])
        average_x = np.mean(array[:, 0])
        average_y = np.mean(array[:, 1])
        mask2center[nucleus] = [average_x, average_y]
        sizes.append(len(mask_dict[nucleus]))

    with open(bin_file, 'r') as f:
        next(f)
        lines = [line.strip().split() for line in f]

    #find xmin ymin
    xmin = float('inf') # 找出gem文件中x坐标最小值
    ymin = float('inf') # 找出gem文件中y坐标最小值
    geneid = {} # 记录每个gene对应的id， {gene:id}
    genecnt = 0 # 记录临时索引
    id2gene = {} #记录每个id对应的gene， {id:gene}

    for line in lines:
        gene, x, y, count = line
        x, y = int(x), int(y)
        xmin = min(xmin, x)
        ymin = min(ymin, y)
        if gene not in geneid:
            geneid[gene] = genecnt
            id2gene[genecnt] = gene
            genecnt += 1

    idx2exp = {} # 每块区域的基因表达量， {idx:{geneid:count}}
    downrs = bin_size
    for line in lines:
        gene, x, y, count = line
        x = int(x) - xmin
        y = int(y) - ymin
        if gene not in geneid:
            continue
        if x < int(startx) or x >= int(startx) + int(patchsizex) or y < int(starty) or y >= int(starty) + int(patchsizey):
            continue
        idx = int(math.floor((x - int(startx)) / downrs) * math.ceil(patchsizey / downrs) + math.floor((int(y) - int(starty)) / downrs))
        if idx not in idx2exp:
            idx2exp[idx] = {}
            idx2exp[idx][geneid[gene]] = int(count)
        elif geneid[gene] not in idx2exp[idx]:
            idx2exp[idx][geneid[gene]] = int(count)
        else:
            idx2exp[idx][geneid[gene]] += int(count)

    del lines

    all_exp_merged_bins = lil_matrix((int(math.ceil(patchsizex / downrs) * math.ceil(patchsizey / downrs)), genecnt), dtype=np.int8)
    for idx in idx2exp:
        for gid in idx2exp[idx]:
            all_exp_merged_bins[idx, gid] = idx2exp[idx][gid]
            #print(idx, gid, idx2exp[idx][gid])
    all_exp_merged_bins = all_exp_merged_bins.tocsr()
    #print(all_exp_merged_bins.shape)

    all_exp_merged_bins_ad = ad.AnnData(
        all_exp_merged_bins,
        obs=pd.DataFrame(index=[i for i in range(all_exp_merged_bins.shape[0])]),
        var=pd.DataFrame(index=[i for i in range(all_exp_merged_bins.shape[1])]),
    )


    sc.pp.highly_variable_genes(all_exp_merged_bins_ad, n_top_genes=2000, flavor='seurat_v3', span=1.0)
    selected_index = all_exp_merged_bins_ad.var[all_exp_merged_bins_ad.var.highly_variable].index
    selected_index = list(selected_index)
    selected_index = [int(i) for i in selected_index]
    #selected_index = geneidx[selected_index]
    #print(selected_index, len(selected_index))

    with open(tmp_dir+'/variable_genes' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.txt', 'w') as fw:
        for id in selected_index:
            fw.write(id2gene[id] + '\n')

    #check total gene counts
    all_exp_merged_bins = all_exp_merged_bins.toarray()[:, selected_index]
    #print(all_exp_merged_bins.shape)
    # cell_bins = {}
    x_train_tmp = []
    x_train = []
    x_train_pos = []
    y_train = []
    y_binary_train = []
    x_train_bg_tmp = []
    x_train_bg = []
    x_train_pos_bg = []
    y_train_bg = []
    y_binary_train_bg = []
    x_test_tmp = []
    x_test= []
    x_test_pos = []
    offsets = []
    for dis in range(1, 11):
        for dy in range(-dis, dis + 1):
            offsets.append([-dis * downrs, dy * downrs])
        for dy in range(-dis, dis + 1):
            offsets.append([dis * downrs, dy * downrs])
        for dx in range(-dis + 1, dis):
            offsets.append([dx * downrs, -dis * downrs])
        for dx in range(-dis + 1, dis):
            offsets.append([dx * downrs, dis * downrs])
    for i in range(adatasub.layers[layer_name].shape[0]):
        if (i + 1) % 100 == 0:
            print("finished {0:.0%}".format(i / adatasub.layers[layer_name].shape[0]))
        for j in range(adatasub.layers[layer_name].shape[1]):
            if (not i % downrs == 0) or (not j % downrs == 0):
                continue
            idx = int(math.floor(i / downrs) * math.ceil(patchsizey / downrs) + math.floor(j / downrs))
            if adatasub.layers[layer_name][i, j] > 0:
                if idx >= 0 and idx < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx, :]) > 0:
                    x_train_sample = [all_exp_merged_bins[idx, :]]
                    x_train_pos_sample = [[i, j]]
                    y_train_sample = [mask2center[adatasub.layers[layer_name][i, j]]]
                    for dx, dy in offsets:
                        if len(x_train_sample) == n_neighbor:
                            break
                        x = i + dx
                        y = j + dy
                        if x < 0 or x >= adatasub.layers[layer_name].shape[0] or y < 0 or y >= adatasub.layers[layer_name].shape[1]:
                            continue
                        idx_nb = int(math.floor(x / downrs) * math.ceil(patchsizey / downrs) + math.floor(y / downrs))
                        if idx_nb >= 0 and idx_nb < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx_nb, :]) > 0:
                            x_train_sample.append(all_exp_merged_bins[idx_nb, :])
                            x_train_pos_sample.append([x, y])
                    if len(x_train_sample) < n_neighbor:
                        continue
                    x_train_tmp.append(x_train_sample)
                    if len(x_train_tmp) > 500:
                        x_train.extend(x_train_tmp)
                        x_train_tmp = []
                        #print(np.array(x_train).shape)
                    #print(x_train)
                    #print(len(x_train_tmp))
                    x_train_pos.append(x_train_pos_sample)
                    y_train.append(y_train_sample)
                    y_binary_train.append(1)
            else:
                if idx >= 0 and idx < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx, :]) > 0:
                    backgroud = True
                    for nucleus in mask2center:
                        if (i - mask2center[nucleus][0]) ** 2 + (j - mask2center[nucleus][1]) ** 2 <= 900 or adatasub.layers['stain'][i, j] > 10:
                            backgroud = False
                            break
                    if backgroud:
                        if len(x_train_bg) + len(x_train_bg_tmp) >= len(x_train) + len(x_train_tmp):
                            continue
                        x_train_sample = [all_exp_merged_bins[idx, :]]
                        x_train_pos_sample = [[i, j]]
                        y_train_sample = [[-1, -1]]
                        for dx, dy in offsets:
                            if len(x_train_sample) == n_neighbor:
                                break
                            x = i + dx
                            y = j + dy
                            if x < 0 or x >= adatasub.layers[layer_name].shape[0] or y < 0 or y >= adatasub.layers[layer_name].shape[1]:
                                continue
                            idx_nb = int(math.floor(x / downrs) * math.ceil(patchsizey / downrs) + math.floor(y / downrs))
                            if idx_nb >= 0 and idx_nb < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx_nb, :]) > 0:
                                x_train_sample.append(all_exp_merged_bins[idx_nb, :])
                                x_train_pos_sample.append([x, y])
                        if len(x_train_sample) < n_neighbor:
                            continue
                        x_train_bg_tmp.append(x_train_sample)
                        if len(x_train_bg_tmp) > 500:
                            x_train_bg.extend(x_train_bg_tmp)
                            x_train_bg_tmp = []
                        #print(len(x_train_bg_tmp))
                        x_train_pos_bg.append(x_train_pos_sample)
                        y_train_bg.append(y_train_sample)
                        y_binary_train_bg.append(0)
                        #print(np.sum(posexp[str(i) + ':' + str(j)]), adatasub.X[i,j])
                    else:
                        x_test_sample = [all_exp_merged_bins[idx, :]]
                        x_test_pos_sample = [[i, j]]
                        for dx, dy in offsets:
                            if len(x_test_sample) == n_neighbor:
                                break
                            x = i + dx
                            y = j + dy
                            exp_merge = np.zeros(len(selected_index))
                            if x < 0 or x >= adatasub.layers[layer_name].shape[0] or y < 0 or y >= adatasub.layers[layer_name].shape[1]:
                                continue
                            idx_nb = int(math.floor(x / downrs) * math.ceil(patchsizey / downrs) + math.floor(y / downrs))
                            if idx_nb >= 0 and idx_nb < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx_nb, :]) > 0:
                                x_test_sample.append(all_exp_merged_bins[idx_nb, :])
                                x_test_pos_sample.append([x, y])
                        if len(x_test_sample) < n_neighbor:
                            continue
                        x_test_tmp.append(x_test_sample)
                        if len(x_test_tmp) > 500:
                            x_test.extend(x_test_tmp)
                            x_test_tmp = []
                        x_test_pos.append(x_test_pos_sample)#
    x_train.extend(x_train_tmp)
    x_train_bg.extend(x_train_bg_tmp)
    x_test.extend(x_test_tmp)

    x_train = np.array(x_train)
    x_train_pos = np.array(x_train_pos)
    y_train = np.vstack(y_train)
    y_binary_train = np.array(y_binary_train)
    x_train_bg = np.array(x_train_bg)
    x_train_pos_bg = np.array(x_train_pos_bg)
    y_train_bg = np.vstack(y_train_bg)
    y_binary_train_bg = np.array(y_binary_train_bg)
    #print(x_train.shape, x_train_pos.shape, y_train.shape, y_binary_train.shape, x_train_bg.shape, x_train_pos_bg.shape, y_train_bg.shape, y_binary_train_bg.shape)

    bg_index = np.arange(len(x_train_bg))
    np.random.shuffle(bg_index)
    x_train = np.vstack((x_train, x_train_bg[bg_index[:len(x_train)]]))
    x_train_pos = np.vstack((x_train_pos, x_train_pos_bg[bg_index[:len(x_train_pos)]]))
    y_train = np.vstack((y_train, y_train_bg[bg_index[:len(y_train)]]))
    y_binary_train = np.hstack((y_binary_train, y_binary_train_bg[bg_index[:len(y_binary_train)]]))
    #print(x_train.shape, x_train_pos.shape, y_train.shape, y_binary_train.shape)

    x_test= np.array(x_test)
    x_test_pos = np.array(x_test_pos)
    #print(x_test.shape, x_test_pos.shape)

    np.savez_compressed(tmp_dir+'/x_train_' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.npz', x_train=x_train)
    np.savez_compressed(tmp_dir+'/x_train_pos_' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.npz', x_train_pos=x_train_pos)
    np.savez_compressed(tmp_dir+'/y_train_' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.npz', y_train=y_train)
    np.savez_compressed(tmp_dir+'/y_binary_train_' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.npz', y_binary_train=y_binary_train)
    np.savez_compressed(tmp_dir+'/x_test_' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.npz', x_test=x_test)
    np.savez_compressed(tmp_dir+'/x_test_pos_' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.npz', x_test_pos=x_test_pos)
