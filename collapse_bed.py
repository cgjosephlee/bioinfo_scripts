#!/usr/bin/env python3

def collapse_spans(total_span, retrun_len=False):
    '''
    0-base coordinates
    make sure end is always larger than start
    '''
    total_span.sort(key=lambda x: (x[0], x[1]))
    final_span = []
    span = total_span[0]
    for i in range(1, len(total_span)):
        # span0 = total_span[i-1]
        span1 = total_span[i]
        if span[1] < span1[0]:
            # not overlap, next
            final_span.append(span)
            span = span1
        elif span[1] < span1[1]:
            # overlap and new end is larger than old end, extend span
            span[1] = span1[1]
    # last one
    final_span.append(span)
    if retrun_len:
        final_span_len = 0
        for x in final_span:
            final_span_len += (x[1] - x[0])
        return final_span_len
    else:
        return final_span

if __name__ == '__main__':
    pass
