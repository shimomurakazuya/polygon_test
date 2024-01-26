/*****************************************************************************/
/**
 *  @file   ArrowGlyph.h
 *  @author Naohisa Sakamoto
 */
/*----------------------------------------------------------------------------
 *
 *  Copyright (c) Visualization Laboratory, Kyoto University.
 *  All rights reserved.
 *  See http://www.viz.media.kyoto-u.ac.jp/kvs/copyright/ for details.
 *
 *  $Id: ArrowGlyph.h 1797 2014-08-04 01:36:37Z naohisa.sakamoto@gmail.com $
 */
/*****************************************************************************/
#ifndef KVS__ARROW_GLYPH_H_INCLUDE
#define KVS__ARROW_GLYPH_H_INCLUDE

#include "GlyphBase.h"
#include <kvs/Module>
#include <kvs/Type>
#include <kvs/OpenGL>
#include <kvs/PolygonObject>
#include <kvs/Quaternion>


namespace kvs
{

class ObjectBase;
class Camera;
class Light;
class VolumeObjectBase;
class TransferFunction;
class RGBColor;

/*===========================================================================*/
/**
 *  @brief  ArrowGlyph class.
 */
/*===========================================================================*/
class ArrowGlyph : public kvs::GlyphBase , public kvs::PolygonObject
{
    kvsModule( kvs::ArrowGlyph, Renderer );
    kvsModuleBaseClass( kvs::GlyphBase );
    kvsModuleSuperClass( kvs::PolygonObject );

public:

    enum ArrowType
    {
        LineArrow = 0,
        TubeArrow
    };

private:

    ArrowType m_type; ///< arrow type
    const kvs::VolumeObjectBase* m_volume; ///< pointer to the volume object (reference)
    kvs::ValueArray<kvs::Real32> m_position;
    //kvs::Vector3f m_position = kvs::Vector3f( 0.f, 0.f, 0.f);
    kvs::Real32   m_size;
    kvs::RGBColor m_color; 
    kvs::UInt8 m_opacity;
    int m_glyph_num;
    kvs::ValueArray<kvs::Real32> m_directions;
    kvs::ValueArray<kvs::Real32> m_coords;
    kvs::ValueArray<kvs::UInt32> m_connection;
    kvs::ValueArray<kvs::Real32> m_normal;


public:

    ArrowGlyph();
    ArrowGlyph( const kvs::VolumeObjectBase* volume );
    ArrowGlyph( const kvs::VolumeObjectBase* volume, const kvs::TransferFunction& transfer_function );
    ArrowGlyph( const int nglyphs, const kvs::ValueArray<kvs::Real32>& position, const kvs::ValueArray<kvs::Real32>& direction);

    void setType( const ArrowType type ) { m_type = type; }
    //void setPosition(const kvs::Vector3f position ){ m_position = position;}
    void setPosition(const kvs::ValueArray<kvs::Real32>& position ){m_position = position;}
    void setSize(const kvs::Real32 size ){ m_size = size; }
    void setRGBColor(const kvs::RGBColor color ){ m_color = color ;}
    void setOpacity(const kvs::UInt8 opacity ){ m_opacity = opacity ;}
    void setGlyphNum(const int glyph_num ) { m_glyph_num = glyph_num;}
    //void setDirection()
    //void ConvertPolygon();
    void setPolygon(); 
    void GenerateNormalizedPolygon();
    void Transform(const int nglyphs,
                   const kvs::ValueArray<kvs::Real32>& position,
                   const kvs::ValueArray<kvs::Real32>& direction,
                   const kvs::Real32 m_size );
    ArrowType type() const { return m_type; }
    void exec( kvs::ObjectBase* object, kvs::Camera* camera, kvs::Light* light );

private:

    void attach_volume( const kvs::VolumeObjectBase* volume );
    void draw();
    void draw_lines();
    void draw_tubes();
    void draw_line_element( const kvs::RGBColor& color, const kvs::UInt8 opacity );
    void draw_tube_element( const kvs::RGBColor& color, const kvs::UInt8 opacity );
    void initialize();
};

} // end of namespace kvs

#endif // KVS__ARROW_GLYPH_H_INCLUDE
