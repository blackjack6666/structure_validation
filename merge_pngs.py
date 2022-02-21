### merge pngs to one horizontally
import sys
from PIL import Image
from glob import glob
from commons import chunks, protein_info_from_fasta

fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
fasta_info_dict = protein_info_from_fasta(fasta_file)

images = glob('D:/data/native_protein_digestion/12072021/low_confident_coverage/*distance.png')
images_slice = [each for each in chunks(images, 7)]
for each in images_slice:
    prt_id = each[0].split('\\')[-1].split('_')[0]
    gene = fasta_info_dict[prt_id][0]
    print(prt_id)
    data_source = each[0].split('\\')[-1].split('_')[-1].split('.png')[0]
    images = [Image.open(x) for x in each]
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]

    new_im.save('D:/data/native_protein_digestion/12072021/low_confident_coverage/' + gene + '_' + data_source + '.png')
