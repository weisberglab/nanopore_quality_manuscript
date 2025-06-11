import pysam
import sys

CHUNK_SIZE = 400

def split_read_into_chunks(read):
    """
    Splits a read into CHUNK_SIZE bp chunks, keeping genomic alignment and adjusting the CIGAR accordingly.
    """
    if read.is_unmapped:
        return []
    #ignore secondary alignments
    if read.flag == 256 or read.flag == 272:
        return []

    seq = read.query_sequence
    qual = read.query_qualities
    cigar = read.cigartuples
    pos = read.reference_start
    # Flatten the query using the CIGAR to get mapping from query pos -> ref pos
    chunks = []
    query_index = 0
    ref_index = pos
    cigar_index = 0
    ops = []
    
    for op, length in cigar:
        if op == 0:  # Match/Mismatch (M)
            ops.extend([('M', ref_index + i, query_index + i) for i in range(length)])
            ref_index += length
            query_index += length
        elif op == 1:  # Insertion (I)
            ops.extend([('I', ref_index, query_index + i) for i in range(length)])
            query_index += length
        elif op == 2:  # Deletion (D)
            ops.extend([('D', ref_index + i, query_index) for i in range(length)])
            ref_index += length
        elif op == 3:  # Skipped region (N)
            ref_index += length
        elif op == 4:  # Soft clip
            query_index += length
        elif op == 5:  # Hard clip
            continue
        elif op == 7 or op == 8:  # = or X
            ops.extend([(op, ref_index + i, query_index + i) for i in range(length)])
            ref_index += length
            query_index += length

    # Group into chunks
    chunk_ops = []
    current_chunk = []
    current_query_len = 0

    for item in ops:
        if current_query_len >= CHUNK_SIZE:
            chunk_ops.append(current_chunk)
            current_chunk = []
            current_query_len = 0
        current_chunk.append(item)
        if item[0] != 'D':  # Only non-deletion operations advance the query
            current_query_len += 1

    if current_chunk:
        chunk_ops.append(current_chunk)

    for i, chunk in enumerate(chunk_ops):
        if not chunk:
            continue

        new_read = read.__copy__()
        new_seq = []
        new_qual = []

        new_cigar = []
        last_op = None
        op_len = 0
        new_start = chunk[0][1]

        for op, ref_pos, q_pos in chunk:
            if op == 'M' or op == 7 or op == 8:
                new_seq.append(seq[q_pos])
                new_qual.append(qual[q_pos])
                op_code = 0 if op == 'M' else op
            elif op == 'I':
                new_seq.append(seq[q_pos])
                new_qual.append(qual[q_pos])
                op_code = 1
            elif op == 'D':
                op_code = 2
            else:
                continue

            if last_op is None or op_code != last_op:
                if last_op is not None:
                    new_cigar.append((last_op, op_len))
                last_op = op_code
                op_len = 1
            else:
                op_len += 1

        if last_op is not None:
            new_cigar.append((last_op, op_len))

        new_read.cigartuples = new_cigar
        new_read.query_sequence = ''.join(new_seq)
        new_read.query_qualities = new_qual
        new_read.reference_start = new_start
        new_read.query_name += f"_chunk{i+1}"
        #print(new_read)
        chunks.append(new_read)

    return chunks

def process_bam(input_bam, output_bam):
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", template=in_bam)

    for read in in_bam:
        if read.is_unmapped:
            continue
        for chunk in split_read_into_chunks(read):
            out_bam.write(chunk)

    in_bam.close()
    out_bam.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python chunk_bam_reads.py input.bam output.bam")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]

    process_bam(input_bam, output_bam)

