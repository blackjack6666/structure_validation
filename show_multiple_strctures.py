import commons
from collections import defaultdict
import re

def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, regex_dict=None):
    """

    :param psm_list:
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """

    from collections import Counter
    from collections import defaultdict
    from heapq import heappush, heappop

    id_freq_array_dict = {}
    ptm_site_counting = defaultdict(int)
    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
    h = []
    regex_pat = '\w{1}\[\d+\.?\d+\]' # universal ptm pattern
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    for each in psm_list:

        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)

    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    ptm_index_line_dict = {each:commons.zero_line_for_seq(seq_line) for each in regex_dict} if regex_dict else False
    separtor_pos_array = commons.separator_pos(seq_line)

    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1]+1]+=len(peptide_psm_dict[matched_pep])
        # ptm assign might need to be optimized
        if ptm_index_line_dict:  # if ptm enabled
            for psm in peptide_psm_dict[matched_pep]:
                for ptm in regex_dict:
                    # only keep one ptm in psm if there are multiple for correct index finding
                    new_psm = re.sub('(?!'+ptm[1:]+')\[\d+\.?\d+\]','',psm)

                    ptm_mod = re.findall(ptm, new_psm)
                    print(ptm_mod)
                    if ptm_mod:
                        for ele in ptm_mod:
                            print (new_psm)
                            ptm_idx = new_psm.find(ele)
                            print(matched_pep, tp[0], ele, ptm_idx)
                            ptm_index_line_dict[ptm][tp[0] + ptm_idx] += 1

    time_start = time.time()
    for i in range(len(separtor_pos_array)-1):
        zero_line_slice = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]]
        percentage_cov = np.count_nonzero(zero_line_slice)/len(zero_line_slice)*100
        heappush(h,(percentage_cov,id_list[i],zero_line_slice))
        id_freq_array_dict[id_list[i]] = zero_line_slice.tolist()


        if ptm_index_line_dict:
            id_ptm_idx_dict[id_list[i]]= {ptm:np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]+1
                                          for ptm in ptm_index_line_dict}
    print (time.time()-time_start)


    return id_freq_array_dict, id_ptm_idx_dict, [heappop(h) for i in range(len(h))][::-1]