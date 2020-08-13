/***************************************************************************
                          bc3d.h  -  description
                             -------------------
    begin                : Mon Nov 17 2003
    copyright            : (C) 2003 by jcch
    email                : jcch@mail.ist.utl.pt
 ***************************************************************************/

// power of 2 always
#define BC_TYPE_NONE                      0x00000000
#define BC_TYPE_FAR_FIELD                 0x00000001

#define BC_TYPE_SUBSONIC_INLET            0x00000010
#define BC_TYPE_SUBSONIC_INLET            0x00000010
#define BC_TYPE_SUBSONIC_INLET_WALL       0x00000020
#define BC_TYPE_SUBSONIC_OUTLET           0x00000040

#define BC_TYPE_SUPERSONIC_INLET          0x00000100
#define BC_TYPE_SUPERSONIC_OUTLET         0x00000400

#define BC_TYPE_WALL                      0x00001000
#define BC_TYPE_MASK                      0x00FFFFFF

#define BC_INDEX_MASK                     0xFF000000
#define BC_INDEX_MASTER                   0x01000000

#define BC_INDEX_SLAVE                    0x02000000
