import pandas as pd
import lsst.geom
import lsst.daf.butler as daf_butler

def get_overlap_region(row, bbox):
    llc = lsst.geom.Point2I(max(row.x0, bbox.minX),
                            max(row.y0, bbox.minY))
    urc = lsst.geom.Point2I(min(row.x0 + row.width - 1, bbox.maxX),
                            min(row.y0 + row.height - 1, bbox.maxY))
    return lsst.geom.Box2I(llc, urc)

def area(row):
    return row.width*row.height

def tabulate_defects(butler, dsref, colthresh=100):
    defects = butler.getDirect(dsref).toDict()
    df0 = pd.DataFrame({_: defects[_] for _ in ('x0', 'y0', 'width', 'height')})
    camera = butler.get('camera', instrument=dsref.dataId['instrument'])
    detector = dsref.dataId['detector']
    det = camera[detector]
    total_area = 0
    for _, row in df0.iterrows():
        total_area += area(row)
    total_region_area = 0
    amp_data = {}
    for i, amp in enumerate(det):
        bbox = amp.getBBox()
        df = df0.query(f'{bbox.minX} - width < x0 <= {bbox.maxX} and '
                       f'{bbox.minY} - height < y0 <= {bbox.maxY}')
        regions = []
        for _, row in df.iterrows():
            regions.append(get_overlap_region(row, bbox))
        total_region_area += sum(area(_) for _ in regions)

        bad_columns = set()
        for region in regions:
            if region.height > colthresh:
                bad_columns.add(region.minX)
        isolated_regions = [_ for _ in regions if _.minX not in bad_columns]

        bad_pixel_area = (sum(area(_) for _ in isolated_regions)
                          + bbox.height*len(bad_columns))
        amp_data[amp.getName()] = (bad_pixel_area, len(bad_columns))
    assert total_area == total_region_area
    return amp_data

if __name__ == '__main__':
    repo = '/repo/main'
    #collections = ['u/jchiang/dark_defects_13162_w_2022_39']
    #colthresh = 100
    collections = ['u/jchiang/bright_defects_13162_w_2022_39']
    colthresh = 20
    butler = daf_butler.Butler(repo, collections=collections)

    dstype = 'defects'
    dsrefs = list(set(butler.registry.queryDatasets(dstype)))
    print("len(dsrefs):", len(dsrefs))

    for dsref in dsrefs[:10]:
        amp_data = tabulate_defects(butler, dsref, colthresh=colthresh)
        print(dsref.dataId['detector'])
        print(amp_data)
        print()
