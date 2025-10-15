import numpy as np
import pandas as pd
import logging
import os
import shutil

class EigenstratReader:

    def __init__(self, filename, chunk_size=1000):
        self.filename = filename
        self.chunk_size = chunk_size
        self.file = None
        self.is_packed = None
        self.n_ind = None
        self.n_snp = None
        self.row_len = None
        self.current_pos = 0
        
    def __iter__(self):
        try:
            with open(self.filename, 'rb') as f:
                header = f.read(100).decode().strip()
                try:
                    header_parts = header.split('\n')[0].split()
                    if len(header_parts) >= 3:
                        self.n_ind = int(header_parts[1])
                        self.n_snp = int(header_parts[2])
                        self.row_len = max(48, int(np.ceil(self.n_ind * 2 / 8)))
                        self.is_packed = True
                        yield from self._read_packed()
                    else:
                        self.is_packed = False
                        yield from self._read_unpacked()
                        
                except (ValueError, IndexError):
                    self.is_packed = False
                    yield from self._read_unpacked()
                    
        except Exception as e:
            logging.error(f"Error reading genotype file {self.filename}: {str(e)}")
            raise
    
    def _read_packed(self):
        try:
            with open(self.filename, 'rb') as f:
                # skip header
                f.seek(self.row_len)

                for i in range(0, self.n_snp, self.chunk_size):
                    current_chunk_size = min(self.chunk_size, self.n_snp - i)
                    bytes_to_read = current_chunk_size * self.row_len

                    chunk_data = f.read(bytes_to_read)
                    if len(chunk_data) == 0:
                        break

                    actual_rows = len(chunk_data) // self.row_len
                    if actual_rows == 0:
                        break

                    chunk = np.frombuffer(chunk_data[:actual_rows * self.row_len], dtype='uint8').reshape(actual_rows, self.row_len)
                    unpacked = np.unpackbits(chunk, axis=1)[:, :(2 * self.n_ind)]
                    genotypes = 2 * unpacked[:, ::2] + unpacked[:, 1::2]
                    genotypes[genotypes == 3] = 9
                    yield genotypes

        except Exception as e:
            logging.error("Error reading packed genotype data")
            raise
    
    def _read_unpacked(self):
        try:
            with open(self.filename, 'r') as f:
                first_line = f.readline().strip()
                n_ind = len(first_line)
                f.seek(0)
                
                chunk_snp_count = 0
                chunk_array = np.zeros((self.chunk_size, n_ind), dtype=np.int8)

                for line in f:
                    line = line.strip()
                    if len(line) == 0:
                        continue

                    line_data = line[:n_ind]
                    genotype_data = np.frombuffer(line_data.encode('ascii'), dtype=np.int8) - 48

                    chunk_array[chunk_snp_count, :len(genotype_data)] = genotype_data
                    chunk_snp_count += 1

                    if chunk_snp_count >= self.chunk_size:
                        yield chunk_array[:chunk_snp_count].copy()
                        chunk_snp_count = 0

                if chunk_snp_count > 0:
                    yield chunk_array[:chunk_snp_count].copy()
            
        except Exception as e:
            logging.error(f"Error reading unpacked genotype data: {e}")
            raise

def read_ind_file(filename):
    ind_data = pd.read_csv(
        filename,
        sep=r'\s+',
        header=None,
        names=['iid', 'sex', 'population'])
    return ind_data

def read_snp_file(filename):
    snp_data = pd.read_csv(
        filename,
        sep=r'\s+',
        header=None,
        names=["snpID", "chr", "cM", "pos", "allele1", "allele2"])
    return snp_data

class OutputHandler:
    def __init__(self, config, ind_data):
        self.config = config
        self.ind_data = ind_data
        self.packed = config.out_type == 'anc'
        self.n_ind = len(ind_data)
        self.total_snps = 0

    def _write_packed_header(self, outname):
        rLen = max(48, int(np.ceil(self.n_ind * 2 / 8)))
        header = np.zeros(rLen, dtype=np.uint8)

        header_str = f"GENO    {self.n_ind}    {0:>10} ".encode()

        if len(header_str) >= rLen:
            raise ValueError(f"Header too large ({len(header_str)} bytes) for buffer ({rLen} bytes)")

        header[:len(header_str)] = np.frombuffer(header_str, dtype=np.uint8)

        with open(outname, 'wb') as f:
            f.write(header.tobytes())

    def _write_packed_geno(self, geno, outname, first_chunk=False):
        try:
            nInd = geno.shape[1]
            nSnp = geno.shape[0]
            self.total_snps += nSnp  
            rLen = max(48, int(np.ceil(nInd * 2 / 8)))

            mode = 'ab'
            with open(outname, mode) as f:
                packed_data = np.zeros((nSnp, rLen), dtype=np.uint8)
                for snp_idx in range(nSnp):
                    genotypes = geno[snp_idx].copy()
                    genotypes[genotypes == 9] = 3
                    bits = np.zeros(2 * nInd, dtype=np.uint8)
                    bits[::2] = genotypes // 2  
                    bits[1::2] = genotypes % 2 
                    
                    packed_bytes = np.packbits(bits)
                    packed_data[snp_idx, :len(packed_bytes)] = packed_bytes
                
                f.write(packed_data.tobytes())
        except Exception as e:
            logging.error(f"Error writing packed genotype data: {e}")
            safe_remove(outname)
            raise

    def finalise_packed_geno(self, outname):
        try:
            rLen = max(48, int(np.ceil(self.n_ind * 2 / 8)))
            header = np.zeros(rLen, dtype=np.uint8)
            header_str = f"GENO    {self.n_ind}    {self.total_snps} ".encode()

            if len(header_str) >= rLen:
                raise ValueError(f"Header text too large for buffer")

            header[:len(header_str)] = np.frombuffer(header_str, dtype=np.uint8)

            with open(outname, 'r+b') as f:
                f.seek(0)
                f.write(header.tobytes())  # overwrite header

            logging.info(f"Finalised packed genotype file with {self.total_snps} SNPs and {self.n_ind} individuals")

        except Exception as e:
            logging.error(f"Error finalising packed genotype file: {e}")

            if os.path.exists(outname):
                try:
                    os.remove(outname)
                except Exception:
                    pass

            raise

    def _write_unpacked(self, chunk, file):
        is_binary_mode = hasattr(file, 'mode') and 'b' in file.mode
        valid_values = np.isin(chunk, [0, 1, 2, 9])

        if not np.all(valid_values):
            invalid_positions = np.where(~valid_values)
            invalid_values = chunk[invalid_positions]
            raise ValueError(f"Invalid genotype values found: {invalid_values}, at positions: {invalid_positions}. Expecting only 0, 1, 2, or 9.")
        
        geno_data = chunk.copy()
        geno_data = np.clip(geno_data, 0, 2)
        geno_data[chunk == 9] = 9  
        ascii_values = np.where(geno_data != 9, geno_data + 48, 57).astype(np.uint8)

        if not is_binary_mode:
            lines = []
            for row in ascii_values:
                line = row.tobytes().decode('ascii') + '\n'
                lines.append(line)
            
            file.write(''.join(lines))
        else:
            result = bytearray()
            for row in ascii_values:
                result.extend(row.tobytes())
                result.extend(b'\n')
            file.write(result)

    def write_geno(self, geno_chunk, first_chunk=True):
        try:
            if self.packed:
                if first_chunk:
                    self._write_packed_header(self.config.out_geno)
                self._write_packed_geno(geno_chunk, self.config.out_geno, first_chunk=False)
            else:
                mode = 'w' if first_chunk else 'a'
                with open(self.config.out_geno, mode, encoding='ascii') as f:
                    self._write_unpacked(geno_chunk, f)
        except Exception as e:
            logging.error(f"Error writing genotype data: {e}")
            safe_remove(self.config.out_geno)
            raise

    def write_snp(self, snp_chunk, outname, first_chunk=True, ref_alt_alleles=None):
        try:
            mode = 'w' if first_chunk else 'a'
            columns = ['snpID', 'chr', 'cM', 'pos', 'allele1', 'allele2']
            snp_chunk.to_csv(outname, sep='\t', index=False, header=False,
                             mode=mode, columns=columns)
        except Exception as e:
            logging.error(f"Error writing SNP file: {e}")
            safe_remove(outname)
            raise

    def write_ind(self, ind_data, output_path):
        try:
            ind_data.to_csv(
                output_path, 
                sep='\t', 
                header=False, 
                index=False,
                columns=['iid', 'sex', 'population'])
        except Exception as e:
            logging.error(f"Error writing IND file: {e}")

            # remove partial file
            safe_remove(output_path)
            raise

def backup_file(path: str, extension: str = ".backup") -> str:
    if not os.path.exists(path):
        logging.warning(f"File for backup does not exist: {path}")
        return None
    backup_path = f"{path}{extension}"
    try:
        shutil.copy2(path, backup_path)
        logging.info(f"Backed up file: {path} to {backup_path}")
        return backup_path
    except Exception as e:
        logging.error(f"Error backing up file {path}: {e}")
        return None

def safe_remove(path: str):
    try:
        if os.path.exists(path):
            os.remove(path)
            logging.info(f"Removed file: {path}")
    except Exception as e:
        logging.warning(f"Error removing file {path}: {e}")