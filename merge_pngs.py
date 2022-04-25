### merge pngs to one horizontally
import sys
from PIL import Image
from glob import glob
from commons import chunks, protein_info_from_fasta

fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta'
fasta_info_dict = protein_info_from_fasta(fasta_file)

images = glob('F:/native_digestion/chymotrypsin_4_16/mapping_3d_trypsin_12072021/*.png')

trypsin_images = glob('F:/native_digestion/chymotrypsin_4_16/mapping_3d_trypsin_12072021/*_combined.png')
chymotrypsin_image = [glob('F:/native_digestion/chymotrypsin_4_16/mapping_3d\\' + each.split('\\')[-1])[0] for each in
                      trypsin_images]

print(len(trypsin_images))
print(len(chymotrypsin_image))

# images_slice = [each for each in chunks(images, 7)]
images_slice = list(zip(trypsin_images, chymotrypsin_image))
print(images_slice)
for each in images_slice:
    print(each)
    # prt_id = each[0].split('\\')[-1].split('_')[0]
    # gene = fasta_info_dict[prt_id][0]
    # print(prt_id)
    # data_source = each[0].split('\\')[-1].split('_')[-1].split('.png')[0]
    images = [Image.open(x) for x in each]
    widths, heights = zip(*(i.size for i in images))

    # total_width = sum(widths)
    # max_height = max(heights)

    ### vertically merge
    total_width = max(widths)
    max_height = sum(heights)

    ### create new blank image
    new_im = Image.new('RGB', (total_width, max_height), color=(255, 255, 255, 0))
    # new_im = Image.new(images[0].mode, (total_width,max_height))

    ### add x offset in horizontal merge mode
    # x_offset = 0
    # for im in images:
    #     new_im.paste(im, (x_offset, 0))
    #     x_offset += im.size[0]

    ### add y offset in vertical merge mode
    y_offset = 0
    for im in images:
        new_im.paste(im, (0, y_offset))
        y_offset += im.height

    new_im.save(
        'F:/native_digestion/chymotrypsin_4_16/mapping_t_ct_combined/' + each[0].split('\\')[-1].split('_')[
            0] + '_t_ct_combined.png')
