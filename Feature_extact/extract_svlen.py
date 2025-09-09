#!/usr/bin/env python3
import argparse
import sys

def parse_info(info_str):
    info = {}
    for item in info_str.split(';'):
        if not item or item == '.':
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            info[k] = v
        else:
            info[item] = True  # flag
    return info

def get_info_value(info_dict, key):
    """大小写不敏感获取 INFO 值；flag -> 'True'；不存在 -> '.'"""
    if key in info_dict:
        v = info_dict[key]
        return 'True' if v is True else str(v)
    low = key.lower()
    for k, v in info_dict.items():
        if k.lower() == low:
            return 'True' if v is True else str(v)
    return '.'

def compute_ac_from_gt(format_str, sample_fields, n_alt):
    """从 GT 计算 AC（按 ALT 顺序，逗号分隔）。"""
    if not format_str or not sample_fields:
        return '.'
    keys = format_str.split(':')
    try:
        gt_idx = keys.index('GT')
    except ValueError:
        return '.'
    counts = [0] * n_alt
    for s in sample_fields:
        if s in ('.', ''):
            continue
        parts = s.split(':')
        if gt_idx >= len(parts):
            continue
        gt = parts[gt_idx]
        if gt in ('.', './.', '.|.'):
            continue
        alleles = gt.replace('|', '/').split('/')
        for a in alleles:
            if a.isdigit():
                ai = int(a)
                if 1 <= ai <= n_alt:
                    counts[ai - 1] += 1
    if any(counts):
        return ','.join(str(c) for c in counts)
    return '0' if n_alt == 1 else ','.join('0' for _ in range(n_alt))

def normalize_extra_keys(extra_args_list):
    """
    将多次 --extra 及逗号分隔统一展开去重，保持顺序。
    输入示例：['vals', 'AF,CIPOS', 'CILEN']
    输出：['vals','AF','CIPOS','CILEN']
    """
    if not extra_args_list:
        return ['vals']  # 默认提取 vals
    seen = set()
    out = []
    for chunk in extra_args_list:
        for k in (x.strip() for x in chunk.split(',') if x.strip()):
            kl = k.lower()
            if kl not in seen:
                seen.add(kl)
                out.append(k)  # 保留原始大小写书写
    return out

def infer_svtype_from_alt(info_dict, alt, sv_type_value):
    """若 INFO 中没给 SVTYPE，尝试从 ALT 判断 BND。"""
    if sv_type_value != '.':
        return sv_type_value
    if '[' in alt or ']' in alt:
        return 'BND'
    return sv_type_value

def compute_svlen(sv_type_value, info, pos):
    """
    长度计算规则：
    - BND: 0
    - 其次用 SVLEN（取首个、取绝对值）
    - 否则用 END-POS+1（不为负）
    - 都没有 -> '.'
    """
    if sv_type_value == 'BND':
        return '0'
    svlen_txt = get_info_value(info, 'SVLEN')
    if svlen_txt != '.':
        svlen_raw = svlen_txt.split(',')[0]
        try:
            return str(abs(int(svlen_raw)))
        except ValueError:
            return '.'
    end_txt = get_info_value(info, 'END')
    if end_txt != '.':
        try:
            end = int(end_txt)
            return str(max(0, end - pos + 1))
        except ValueError:
            return '.'
    return '.'

def extract_sv_info(vcf_file, output_file, extra_keys):
    with open(vcf_file, 'r') as f, open(output_file, 'w') as out:
        # 写表头
        header_cols = ['CHROM', 'POS', 'SVLEN', 'SVTYPE', 'AC'] + extra_keys
        out.write('\t'.join(header_cols) + '\n')

        for line in f:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue

            ref = fields[3]
            alt = fields[4]
            info = parse_info(fields[7])

            # SVTYPE：优先 INFO，缺失则从 ALT 判断 BND
            sv_type_value = get_info_value(info, 'SVTYPE')
            sv_type_value = infer_svtype_from_alt(info, alt, sv_type_value)

            # SVLEN
            sv_len_value = compute_svlen(sv_type_value, info, pos)

            # AC：优先 INFO.AC，否则从 GT 计算
            ac_value = get_info_value(info, 'AC')
            if ac_value == '.':
                n_alt = 0 if alt == '.' else len(alt.split(','))
                format_str = fields[8] if len(fields) > 8 else ''
                sample_fields = fields[9:] if len(fields) > 9 else []
                ac_value = compute_ac_from_gt(format_str, sample_fields, n_alt) if n_alt > 0 else '.'

            # 额外 INFO
            extras = [get_info_value(info, k) for k in extra_keys]

            out.write('\t'.join(map(str, [chrom, pos, sv_len_value, sv_type_value, ac_value] + extras)) + '\n')

    print(f"Output written to {output_file}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Extract SV position, length, type, AC, and extra INFO keys from a VCF."
    )
    parser.add_argument("vcf_file", type=str, help="Path to the input VCF file")
    parser.add_argument("output_file", type=str, help="Path to the output file")
    # 新：--extra 可多次传入，也支持逗号分隔；默认 ['vals']
    parser.add_argument(
        "-E", "--extra",
        action="append",
        default=None,
        help="Extra INFO keys to extract. Can be repeated and/or comma-separated, e.g. "
             "`-E vals -E AF,CIPOS -E CILEN`. Case-insensitive."
    )
    # 兼容旧参数名（如果你脚本里已经写了 --extra-info-keys）
    parser.add_argument(
        "--extra-info-keys",
        type=str,
        default=None,
        help="(Deprecated) Comma-separated INFO keys to extract additionally. "
             "Use -E/--extra instead."
    )
    args = parser.parse_args()

    # 兼容：把旧参数并入新参数
    extra_args = args.extra[:] if args.extra else []
    if args.extra_info_keys:
        extra_args.append(args.extra_info_keys)

    extra_keys = normalize_extra_keys(extra_args)
    extract_sv_info(args.vcf_file, args.output_file, extra_keys)

if __name__ == "__main__":
    main()
