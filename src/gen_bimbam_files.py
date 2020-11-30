# from https://stackoverflow.com/a/312464
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


if __name__ == "__main__":
    import os

    res_directory = "/home/nickr/fastq"

    IDS_NUM = 10000
    SNPS_NUM = 20000
    BLOCK_SIZE = 5000  # There will be one chromosome for 5_000 SNPS

    def create_path(filename): return os.path.join(res_directory, filename)

    snps_pos_file = open(create_path("test_data_snps.txt"), 'w')
    i = 0
    cur_chr = 0
    for snps_on_one_chromosome_range in chunks(range(0, SNPS_NUM), BLOCK_SIZE):
        for _ in snps_on_one_chromosome_range:
            snp_pos_record = '{0}\t{0}\t{1}\n'.format(i, cur_chr)
            snps_pos_file.write(snp_pos_record)
            i += 1
        cur_chr += 1
    import random

   "".join(str(random.uniform(0, 2)+ ", ") for _ in range(N))
    id_str = (str(random.uniform(0, 2)) + ", ") * (IDS_NUM - 1) + '1\n'
    geno_file = open(create_path("test_data_geno.txt"), 'w')
    for snp_i in range(0, SNPS_NUM):
        geno_file.write('{0} X, Y, '.format(snp_i) + id_str)
