import pandas as pd
import lsst.geom
import lsst.daf.butler as daf_butler


__all__ = ['tabulate_defects']


def get_overlap_region(row, bbox):
    llc = lsst.geom.Point2I(max(row.x0, bbox.minX),
                            max(row.y0, bbox.minY))
    urc = lsst.geom.Point2I(min(row.x0 + row.width - 1, bbox.maxX),
                            min(row.y0 + row.height - 1, bbox.maxY))
    return lsst.geom.Box2I(llc, urc)


def tabulate_defects(butler, dsref, colthresh=100):
    defects = butler.getDirect(dsref).toDict()
    df0 = pd.DataFrame({_: defects[_] for _ in ('x0', 'y0', 'width', 'height')})
    camera = butler.get('camera', instrument=dsref.dataId['instrument'])
    detector = dsref.dataId['detector']
    det = camera[detector]
    total_area = 0
    for _, row in df0.iterrows():
        total_area += row.width*row.height
    total_region_area = 0
    bad_columns = {}
    bad_pixels = {}
    for i, amp in enumerate(det):
        bbox = amp.getBBox()
        df = df0.query(f'{bbox.minX} - width < x0 <= {bbox.maxX} and '
                       f'{bbox.minY} - height < y0 <= {bbox.maxY}')
        regions = []
        for _, row in df.iterrows():
            regions.append(get_overlap_region(row, bbox))
        total_region_area += sum(_.area for _ in regions)

        bad_cols = set()
        for region in regions:
            if region.height > colthresh:
                bad_cols.add(region.minX)
        bad_columns[amp.getName()] = len(bad_cols)

        isolated_regions = [_ for _ in regions if _.minX not in bad_cols]
        bad_pixel_area = (sum(_.area for _ in isolated_regions)
                          + bbox.height*len(bad_cols))
        bad_pixels[amp.getName()] = bad_pixel_area
    assert total_area == total_region_area
    return bad_columns, bad_pixels
