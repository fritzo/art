#!/usr/bin/env python

import svgwrite
from parsable import parsable


@parsable
def mandala(filename='poirot.svg'):
    '''
    Draw a mandala, inspired by Danna's.
    '''
    d = svgwrite.Drawing(filename)
    d.add(d.line((0, 0), (10, 0), stroke=svgwrite.rgb(10, 10, 16, '%')))
    d.save()


if __name__ == '__main__':
    parsable()
