/*****************************************************************************/
/**
 *  @file   ArrowGlyph.cpp
 *  @author Naohisa Sakamoto
 */
/*----------------------------------------------------------------------------
 *
 *  Copyright (c) Visualization Laboratory, Kyoto University.
 *  All rights reserved.
 *  See http://www.viz.media.kyoto-u.ac.jp/kvs/copyright/ for details.
 *
 *  $Id: ArrowGlyph.cpp 1797 2014-08-04 01:36:37Z naohisa.sakamoto@gmail.com $
 */
/*****************************************************************************/
#include "ArrowGlyph.h"
#include <kvs/OpenGL>
#include <kvs/Quaternion>
#include <kvs/IgnoreUnusedVariable>

namespace
{

const kvs::Real32 LineVertices[12] =
{
     0.0f, 1.0f, 0.0f,
     0.0f, 0.0f, 0.0f,
    -0.2f, 0.8f, 0.0f,
     0.2f, 0.8f, 0.0f
};
const kvs::UInt32 LineConnections[6] =
{
    0, 1,
    0, 2,
    0, 3
};

const kvs::Vec3 ConeTranslation = kvs::Vec3( 0.0f, 0.0f, 0.7f );
inline void DrawCone( const kvs::Vec3& t )
{
    const GLdouble base = 0.15;
    const GLdouble top = 0.0;
    const GLdouble height = 0.3;
    const GLint slices = 20;
    const GLint stacks = 5;
    kvs::OpenGL::Translate( t );
    kvs::OpenGL::DrawCylinder( base, top, height, slices, stacks );
}

const kvs::Vec3 CylinderTranslation = kvs::Vec3( 0.0f, 0.0f, 0.0f );
inline void DrawCylinder( const kvs::Vec3& t )
{
    const GLdouble base = 0.07;
    const GLdouble top = 0.07;
    const GLdouble height = 0.7;
    const GLint slices = 20;
    const GLint stacks = 2;
    kvs::OpenGL::Translate( t );
    kvs::OpenGL::DrawCylinder( base, top, height, slices, stacks );
}

}; // end of namespace


namespace kvs
{

/*===========================================================================*/
/**
 *  @brief  Constructs a new ArrowGlyph class.
 */
/*===========================================================================*/
ArrowGlyph::ArrowGlyph():
    kvs::GlyphBase(),
    m_type( ArrowGlyph::LineArrow ),
    m_volume( NULL )
{
}

/*===========================================================================*/
/**
 *  @brief  Constructs a new ArrowGlyph class.
 *  @param  volume [in] pointer to the volume object
 */
/*===========================================================================*/
ArrowGlyph::ArrowGlyph( const kvs::VolumeObjectBase* volume ):
    kvs::GlyphBase(),
    m_type( ArrowGlyph::LineArrow ),
    m_volume( NULL )
{
    this->attach_volume( volume );
}

/*===========================================================================*/
/**
 *  @brief  Constructs a new ArrowGlyph class.
 *  @param  volume [in] pointer to the Volume object
 *  @param  transfer_function [in] transfer function
 */
/*===========================================================================*/
ArrowGlyph::ArrowGlyph(
    const kvs::VolumeObjectBase* volume,
    const kvs::TransferFunction& transfer_function ):
    kvs::GlyphBase(),
    m_type( ArrowGlyph::LineArrow )
{
    BaseClass::setTransferFunction( transfer_function );
    this->attach_volume( volume );
}

/*===========================================================================*/
/**
 *  @brief  Constructs a new SphereGlyph class.
 *  @param  volume [in] number of glyphs & positions
 */
/*===========================================================================*/
ArrowGlyph::ArrowGlyph( const int nglyphs, const kvs::ValueArray<kvs::Real32>& position, const kvs::ValueArray<kvs::Real32>& direction ):
    kvs::GlyphBase(),
    m_size(1.0f),
    m_opacity( 255),
    m_color(kvs::RGBColor(255, 0, 0 )),
    m_glyph_num(nglyphs)
{
    this->setType( kvs::ArrowGlyph::TubeArrow );
    this->setNormalType( kvs::ArrowGlyph::VertexNormal );
    
    this->GenerateNormalizedPolygon();
    this->Transform(nglyphs, position, direction, m_color, m_size ); //  rotate_scaling_tanslate;
}


/*===========================================================================*/
/**
 *  @brief  Executes rendering process.
 *  @param  object [in] pointer to the volume object
 *  @param  camera [in] pointer to the camera
 *  @param  light [in] pointer to the light
 */
/*===========================================================================*/
void ArrowGlyph::exec( kvs::ObjectBase* object, kvs::Camera* camera, kvs::Light* light )
{
    kvs::IgnoreUnusedVariable( light );
    kvs::IgnoreUnusedVariable( camera );

    const kvs::VolumeObjectBase* volume = kvs::VolumeObjectBase::DownCast( object );
    if ( !volume ) { kvsMessageError("Input object is not volume dat."); return; }
    if ( m_volume != volume ) { this->attach_volume( volume ); }

    BaseClass::startTimer();

    kvs::OpenGL::WithPushedAttrib attrib( GL_CURRENT_BIT | GL_ENABLE_BIT );
    kvs::OpenGL::Enable( GL_DEPTH_TEST );
    this->initialize();
    this->draw();

    BaseClass::stopTimer();
}

/*===========================================================================*/
/**
 *  @brief  Attaches a volume object.
 *  @param  volume [in] pointer to the volume object
 */
/*===========================================================================*/
void ArrowGlyph::attach_volume( const kvs::VolumeObjectBase* volume )
{
    m_volume = volume;

    BaseClass::calculateCoords( volume );

    const std::type_info& type = volume->values().typeInfo()->type();
    if ( type == typeid( kvs::Int8 ) )
    {
        BaseClass::calculateSizes<kvs::Int8>( volume );
        BaseClass::calculateDirections<kvs::Int8>( volume );
        BaseClass::calculateColors<kvs::Int8>( volume );
        BaseClass::calculateOpacities<kvs::Int8>( volume );
    }
    else if ( type == typeid( kvs::Int16 ) )
    {
        BaseClass::calculateSizes<kvs::Int16>( volume );
        BaseClass::calculateDirections<kvs::Int16>( volume );
        BaseClass::calculateColors<kvs::Int16>( volume );
        BaseClass::calculateOpacities<kvs::Int16>( volume );
    }
    else if ( type == typeid( kvs::Int32 ) )
    {
        BaseClass::calculateSizes<kvs::Int32>( volume );
        BaseClass::calculateDirections<kvs::Int32>( volume );
        BaseClass::calculateColors<kvs::Int32>( volume );
        BaseClass::calculateOpacities<kvs::Int32>( volume );
    }
    else if ( type == typeid( kvs::Int64 ) )
    {
        BaseClass::calculateSizes<kvs::Int64>( volume );
        BaseClass::calculateDirections<kvs::Int64>( volume );
        BaseClass::calculateColors<kvs::Int64>( volume );
        BaseClass::calculateOpacities<kvs::Int64>( volume );
    }
    else if ( type == typeid( kvs::UInt8  ) )
    {
        BaseClass::calculateSizes<kvs::UInt8>( volume );
        BaseClass::calculateDirections<kvs::UInt8>( volume );
        BaseClass::calculateColors<kvs::UInt8>( volume );
        BaseClass::calculateOpacities<kvs::UInt8>( volume );
    }
    else if ( type == typeid( kvs::UInt16 ) )
    {
        BaseClass::calculateSizes<kvs::UInt16>( volume );
        BaseClass::calculateDirections<kvs::UInt16>( volume );
        BaseClass::calculateColors<kvs::UInt16>( volume );
        BaseClass::calculateOpacities<kvs::UInt16>( volume );
    }
    else if ( type == typeid( kvs::UInt32 ) )
    {
        BaseClass::calculateSizes<kvs::UInt32>( volume );
        BaseClass::calculateDirections<kvs::UInt32>( volume );
        BaseClass::calculateColors<kvs::UInt32>( volume );
        BaseClass::calculateOpacities<kvs::UInt32>( volume );
    }
    else if ( type == typeid( kvs::UInt64 ) )
    {
        BaseClass::calculateSizes<kvs::UInt64>( volume );
        BaseClass::calculateDirections<kvs::UInt64>( volume );
        BaseClass::calculateColors<kvs::UInt64>( volume );
        BaseClass::calculateOpacities<kvs::UInt64>( volume );
    }
    else if ( type == typeid( kvs::Real32 ) )
    {
        BaseClass::calculateSizes<kvs::Real32>( volume );
        BaseClass::calculateDirections<kvs::Real32>( volume );
        BaseClass::calculateColors<kvs::Real32>( volume );
        BaseClass::calculateOpacities<kvs::Real32>( volume );
    }
    else if ( type == typeid( kvs::Real64 ) )
    {
        BaseClass::calculateSizes<kvs::Real64>( volume );
        BaseClass::calculateDirections<kvs::Real64>( volume );
        BaseClass::calculateColors<kvs::Real64>( volume );
        BaseClass::calculateOpacities<kvs::Real64>( volume );
    }
}

/*===========================================================================*/
/**
 *  @brief  Draw the arrow glyph.
 */
/*===========================================================================*/
void ArrowGlyph::draw()
{
    switch ( m_type )
    {
    case LineArrow: this->draw_lines(); break;
    case TubeArrow: this->draw_tubes(); break;
    default: break;
    }
}

/*===========================================================================*/
/**
 *  @brief  Draw the arrow glyph as lines.
 */
/*===========================================================================*/
void ArrowGlyph::draw_lines()
{
    const size_t npoints = BaseClass::coords().size() / 3;
    const kvs::ValueArray<kvs::Real32> coords = BaseClass::coords();
    const kvs::ValueArray<kvs::UInt8> colors = BaseClass::colors();
    const kvs::ValueArray<kvs::Real32> sizes = BaseClass::sizes();
    const kvs::ValueArray<kvs::UInt8> opacities = BaseClass::opacities();

    if ( BaseClass::directions().size() == 0 )
    {
        for ( size_t i = 0, index = 0; i < npoints; i++, index += 3 )
        {
            const kvs::Vector3f position( coords.data() + index );
            const kvs::Real32 size = sizes[i];
            const kvs::RGBColor color( colors.data() + index );
            const kvs::UInt8 opacity = opacities[i];
            kvs::OpenGL::PushMatrix();
            {
                BaseClass::transform( position, size );
                this->draw_line_element( color, opacity );
            }
            kvs::OpenGL::PopMatrix();
        }
    }
    else
    {
        for( size_t i = 0, index = 0; i < npoints; i++, index += 3 )
        {
            const kvs::Vector3f position( coords.data() + index );
            const kvs::Vector3f direction( BaseClass::directions().data() + index );
            const kvs::Real32 size = sizes[i];
            const kvs::RGBColor color( colors.data() + index );
            const kvs::UInt8 opacity = opacities[i];
            if ( direction.length() > 0.0f )
            {
                kvs::OpenGL::PushMatrix();
                {
                    BaseClass::transform( position, direction, size );
                    this->draw_line_element( color, opacity );
                }
                kvs::OpenGL::PopMatrix();
            }
        }
    }
}

/*===========================================================================*/
/**
 *  @brief  Draw the arrow glyph as polygons.
 */
/*===========================================================================*/
void ArrowGlyph::draw_tubes()
{
    const size_t npoints = BaseClass::coords().size() / 3;
    const kvs::ValueArray<kvs::Real32> coords = BaseClass::coords();
    const kvs::ValueArray<kvs::UInt8> colors = BaseClass::colors();
    const kvs::ValueArray<kvs::Real32> sizes = BaseClass::sizes();
    const kvs::ValueArray<kvs::UInt8> opacities = BaseClass::opacities();

    if ( BaseClass::directions().size() == 0 )
    {
//        for ( size_t i = 0, index = 0; i < npoints; i++, index += 3 )
//        {
//            const kvs::Vector3f position( coords.data() + index );
//            const kvs::Real32 size = sizes[i];
//            const kvs::RGBColor color( colors.data() + index );
//            const kvs::UInt8 opacity = opacities[i];
//            kvs::OpenGL::PushMatrix();
//            {
//                BaseClass::transform( position, size );
//                this->draw_tube_element( color, opacity );
//            }
//            kvs::OpenGL::PopMatrix();
//        }
        for ( size_t i = 0, index = 0; i < npoints; i++, index += 3 )
        {
            const kvs::Vector3f position( coords.data() + index );
            const kvs::Real32 size = sizes[i];
            const kvs::RGBColor color( colors.data() + index );
            const kvs::UInt8 opacity = opacities[i];
            kvs::OpenGL::PushMatrix();
            {
                BaseClass::transform( position, size );
                this->draw_tube_element( color, opacity );
            }
            kvs::OpenGL::PopMatrix();
        }
    }
    else
    {
        for( size_t i = 0, index = 0; i < npoints; i++, index += 3 )
        {
            const kvs::Vector3f position( coords.data() + index );
            const kvs::Vector3f direction( BaseClass::directions().data() + index );
            const kvs::Real32 size = sizes[i];
            const kvs::RGBColor color( colors.data() + index );
            const kvs::UInt8 opacity = opacities[i];
            if ( direction.length() > 0.0f )
            {
                kvs::OpenGL::PushMatrix();
                {
                    BaseClass::transform( position, direction, size );
                    this->draw_tube_element( color, opacity );
                }
                kvs::OpenGL::PopMatrix();
            }
        }
    }
}

/*===========================================================================*/
/**
 *  @brief  Draw the line element.
 *  @param  color [in] color value
 *  @param  opacity [in] opacity value
 */
/*===========================================================================*/
void ArrowGlyph::draw_line_element( const kvs::RGBColor& color, const kvs::UInt8 opacity )
{
    kvs::OpenGL::Begin( GL_LINES );
    kvs::OpenGL::Color( color.r(), color.g(), color.b(), opacity );
    for ( size_t i = 0; i < 6; i++ )
    {
        kvs::OpenGL::Vertex3( ::LineVertices + ::LineConnections[i] * 3 );
    }
    kvs::OpenGL::End();
}

/*===========================================================================*/
/**
 *  @brief  Draw the tube element.
 *  @param  color [in] color value
 *  @param  opacity [in] opacity value
 */
/*===========================================================================*/
void ArrowGlyph::draw_tube_element( const kvs::RGBColor& color, const kvs::UInt8 opacity )
{
    const kvs::Real32 R = -90.0f; // rotation angle
    const kvs::Vec3 V( 1.0f, 0.0f, 0.0f ); // rotation vector
    const kvs::Vec3 T0( 0.0f, 0.0f, 0.7f ); // translation vector (cone)
    const kvs::Vec3 T1( 0.0f, 0.0f, 0.0f ); // translation vector (cylinder)

    kvs::OpenGL::PushMatrix();
    kvs::OpenGL::Rotate( R, V );
    {
        kvs::OpenGL::Color( color.r(), color.g(), color.b(), opacity );

        // Cone.
        kvs::OpenGL::PushMatrix();
        ::DrawCone( T0 );
        kvs::OpenGL::PopMatrix();

        // Cylinder.
        kvs::OpenGL::PushMatrix();
        ::DrawCylinder( T1 );
        kvs::OpenGL::PopMatrix();
    }
    kvs::OpenGL::PopMatrix();
}


void  ArrowGlyph::GenerateNormalizedPolygon()
{

    const kvs::Vec3 T0( 0.0f, 0.0f, 0.7f ); // translation vector (cone)
                                            //    const kvs::Vec3 T1( 0.0f, 0.0f, 0.0f ); // translation vector (cylinder)
    int  n_coords_cone= 126;
    int  n_coords_tube= 63; 
    int  n_coords     = 189; // =(20+1)*(5+1) + (20+1)*(2+1)
    int  n_cells      = 294; // =2*((20+1)*5+(20+1)*2)
    int  n_connection = 3 * n_cells;          //tri cell type

    kvs::ValueArray<kvs::Real32> local_coords;
    kvs::ValueArray<kvs::UInt32> local_connection;
    kvs::ValueArray<kvs::Real32> local_normal;

    local_coords.allocate(3*n_coords);
    local_normal.allocate(3*n_coords);
    local_connection.allocate(n_connection);
    m_coords.allocate(3*n_coords);
    m_normal.allocate(3*n_coords);
    m_connection.allocate(n_connection); 

    int coord_id = 0;
    int connection_id = 0;

    // Cone.
    //        ::DrawCone( T0 );
    //    kvs::OpenGL::Translate( t ); //use T0(( 0.0f, 0.0f, 0.7f ))
    //    kvs::OpenGL::DrawCylinder( base, top, height, slices, stacks );
    GLdouble base = 0.15;
    GLdouble top = 0.0;
    GLdouble height = 0.3;
    GLint slices = 20;
    GLint stacks = 5;

    const int CacheSize = 240;
    const float Pi = 3.14159265358979323846;

    GLint i,j;
    GLfloat sinCache[CacheSize];
    GLfloat cosCache[CacheSize];
    GLfloat sinCache2[CacheSize];
    GLfloat cosCache2[CacheSize];
    GLfloat sinCache3[CacheSize];
    GLfloat cosCache3[CacheSize];
    GLfloat angle;
    GLfloat zLow, zHigh;
    GLfloat length;
    GLfloat deltaRadius;
    GLfloat zNormal;
    GLfloat xyNormalRatio;
    GLfloat radiusLow, radiusHigh;
    int needCache2, needCache3;

    if (slices >= CacheSize) slices = CacheSize-1;

    if (slices < 2 || stacks < 1 || base < 0.0 || top < 0.0 || height < 0.0)
    {
        kvsMessageError("Invalid value.");
        return;
    }

    /* Compute length (needed for normal calculations) */
    deltaRadius = base - top;
    length = std::sqrt(deltaRadius*deltaRadius + height*height);
    if ( length == 0.0 )
    {
        kvsMessageError("Invalid value.");
        return;
    }

    /* Cache is the vertex locations cache */
    /* Cache2 is the various normals at the vertices themselves */
    /* Cache3 is the various normals for the faces */
    needCache2 = 1;
    needCache3 = 0;

    zNormal = deltaRadius / length;
    xyNormalRatio = height / length;

    for (i = 0; i < slices; i++)
    {
        angle = 2 * Pi * i / slices;
        if (needCache2)
        {
            sinCache2[i] = xyNormalRatio * std::sin(angle);
            cosCache2[i] = xyNormalRatio * std::cos(angle);
        }
        sinCache[i] = std::sin(angle);
        cosCache[i] = std::cos(angle);
    }

    if (needCache3)
    {
        for (i = 0; i < slices; i++)
        {
            angle = 2 * Pi * (i-0.5) / slices;
            sinCache3[i] = xyNormalRatio * std::sin(angle);
            cosCache3[i] = xyNormalRatio * std::cos(angle);
        }
    }

    sinCache[slices] = sinCache[0];
    cosCache[slices] = cosCache[0];
    if (needCache2)
    {
        sinCache2[slices] = sinCache2[0];
        cosCache2[slices] = cosCache2[0];
    }
    if (needCache3)
    {
        sinCache3[slices] = sinCache3[0];
        cosCache3[slices] = cosCache3[0];
    }

    /* Note:
     ** An argument could be made for using a TRIANGLE_FAN for the end
     ** of the cylinder of either radii is 0.0 (a cone).  However, a
     ** TRIANGLE_FAN would not work in smooth shading mode (the common
     ** case) because the normal for the apex is different for every
     ** triangle (and TRIANGLE_FAN doesn't let me respecify that normal).
     ** Now, my choice is GL_TRIANGLES, or leave the GL_QUAD_STRIP and
     ** just let the GL trivially reject one of the two triangles of the
     ** QUAD.  GL_QUAD_STRIP is probably faster, so I will leave this code
     ** alone.
     */
    for (j = 0; j <= stacks; j++)
    {
        zLow = j * height / stacks;
        radiusLow = base - deltaRadius * ((float) j / stacks);

        for (i = 0; i <= slices; i++)
        {
            local_coords.at(coord_id   ) =  T0.x() + radiusLow  * sinCache[i];
            local_coords.at(coord_id +1) =  T0.y() + radiusLow  * cosCache[i];
            local_coords.at(coord_id +2) =  T0.z() + zLow;

            local_normal.at(coord_id   ) = sinCache2[i] ;
            local_normal.at(coord_id +1) = cosCache2[i] ;
            local_normal.at(coord_id +2) = zNormal;
            coord_id      +=3;
        }
    }
    int line_size = slices + 1;
    int vertex_index = 0;  
    for (j = 0; j < stacks; j++)
    {
        for (i = 0; i <= slices; i++)
        {
            const int local_vertex_index[4] =
            {   
                vertex_index,
                vertex_index + 1,
                vertex_index + line_size,
                vertex_index + line_size + 1
            }; 
            local_connection.at(connection_id  ) = local_vertex_index[0]; 
            local_connection.at(connection_id+1) = local_vertex_index[2]; 
            local_connection.at(connection_id+2) = local_vertex_index[1]; 

            local_connection.at(connection_id+3) = local_vertex_index[2]; 
            local_connection.at(connection_id+4) = local_vertex_index[3]; 
            local_connection.at(connection_id+5) = local_vertex_index[1]; 
            connection_id += 6;
            vertex_index++;  
        }
    }

    // Cylinder.
    //::DrawCylinder( T1 ) // use T1( 0.0f, 0.0f, 0.0f );
    base = 0.07;
    top = 0.07;
    height = 0.7;
    slices = 20;
    stacks = 2;

    /* Compute length (needed for normal calculations) */
    deltaRadius = base - top;
    length = std::sqrt(deltaRadius*deltaRadius + height*height);
    if ( length == 0.0 )
    {
        kvsMessageError("Invalid value.");
        return;
    }

    /* Cache is the vertex locations cache */
    /* Cache2 is the various normals at the vertices themselves */
    /* Cache3 is the various normals for the faces */
    needCache2 = 1;
    needCache3 = 0;

    zNormal = deltaRadius / length;
    xyNormalRatio = height / length;

    for (i = 0; i < slices; i++)
    {
        angle = 2 * Pi * i / slices;
        if (needCache2)
        {
            sinCache2[i] = xyNormalRatio * std::sin(angle);
            cosCache2[i] = xyNormalRatio * std::cos(angle);
        }
        sinCache[i] = std::sin(angle);
        cosCache[i] = std::cos(angle);
    }

    if (needCache3)
    {
        for (i = 0; i < slices; i++)
        {
            angle = 2 * Pi * (i-0.5) / slices;
            sinCache3[i] = xyNormalRatio * std::sin(angle);
            cosCache3[i] = xyNormalRatio * std::cos(angle);
        }
    }

    sinCache[slices] = sinCache[0];
    cosCache[slices] = cosCache[0];
    if (needCache2)
    {
        sinCache2[slices] = sinCache2[0];
        cosCache2[slices] = cosCache2[0];
    }
    if (needCache3)
    {
        sinCache3[slices] = sinCache3[0];
        cosCache3[slices] = cosCache3[0];
    }

    for (j = 0; j <= stacks; j++)
    {
        zLow = j * height / stacks;
        radiusLow = base - deltaRadius * ((float) j / stacks);

        for (i = 0; i <= slices; i++)
        {
            local_coords.at(coord_id   ) =  radiusLow  * sinCache[i];
            local_coords.at(coord_id +1) =  radiusLow  * cosCache[i];
            local_coords.at(coord_id +2) =  zLow;

            local_normal.at(coord_id   ) = sinCache2[i] ;
            local_normal.at(coord_id +1) = cosCache2[i] ;
            local_normal.at(coord_id +2) = zNormal;
            coord_id      +=3;
        }
    }
    //        line_size = slices + 1;
    vertex_index = 126; 
    for (j = 0; j < stacks; j++)
    {
        for (i = 0; i <= slices; i++)
        {
            const int local_vertex_index[4] =
            {   
                vertex_index,
                vertex_index + 1,
                vertex_index + line_size,
                vertex_index + line_size + 1
            }; 
            local_connection.at(connection_id  ) = local_vertex_index[0]; 
            local_connection.at(connection_id+1) = local_vertex_index[2]; 
            local_connection.at(connection_id+2) = local_vertex_index[1]; 

            local_connection.at(connection_id+3) = local_vertex_index[2]; 
            local_connection.at(connection_id+4) = local_vertex_index[3]; 
            local_connection.at(connection_id+5) = local_vertex_index[1]; 
            connection_id += 6;
            vertex_index++;  
        }
    }

    m_coords = local_coords;
    m_normal = local_normal;
    m_connection = local_connection; 
#ifdef DEBUG
    std::cout << "connection_id = " << connection_id <<std::endl;
    std::cout << "coord_id = " << coord_id <<std::endl;
    std::cout << "object->coord()  = " << float(local_coords[0]) << ", " <<  float(local_coords[1]) << ", " << float(local_coords[2]) << ", "  << 
        float(local_coords[402]) << ", " <<  float(local_coords[403]) << ", " << float(local_coords[404]) << ", " <<
        float(local_coords[564]) << ", " <<  float(local_coords[565]) << ", " << float(local_coords[566]) << ", " << std::endl;
    std::cout << "connection[id]  = " << float(local_connection[0]) << ", " <<  float(local_connection[1]) << ", " << float(local_connection[2]) << ", "  << 
        float(local_connection[879]) << ", " <<  float(local_connection[880]) << ", " << float(local_connection[881]) << ", " << std::endl;
//
//    SuperClass::setCoords( m_coords  ); 
//    SuperClass::setConnections( m_connection );
//    SuperClass::setColor( m_color );
//    SuperClass::setNormals( m_normal );
//    SuperClass::setOpacity( int(m_opacity) );
//    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
//    SuperClass::setColorType( kvs::PolygonObject::PolygonColor );
#endif   
}
/*===========================================================================*/
/**
 *  @brief  translate & scaleing & rotate  arrow glyph.
 */
/*===========================================================================*/
void ArrowGlyph::Transform(const int nglyphs,
                           const kvs::ValueArray<kvs::Real32>& position,
                           const kvs::ValueArray<kvs::Real32>& direction,
                           const kvs::RGBColor color,
                           const kvs::Real32 m_size )
{

    const int n_coords = m_coords.size();
    const int n_connection = m_connection.size();

    kvs::ValueArray<kvs::Real32> gl_coords;
    kvs::ValueArray<kvs::UInt32> gl_connection;
    kvs::ValueArray<kvs::Real32> gl_normal;
    gl_coords.allocate(nglyphs*n_coords);
    gl_normal.allocate(nglyphs*n_coords);
    gl_connection.allocate(nglyphs*n_connection); 

   
    const float Pi = 3.14159265358979323846;
    const kvs::Vec3 DefaultDirection = kvs::Vec3( 0.0, 1.0, 0.0 );
    const kvs::Real32 R = -0.5*Pi; //  -90.0f;              // rotation angle
    const kvs::Vector3f V( 1.0f, 0.0f, 0.0f ); // rotation vector

//    float R[2]={0.f, float(-0.5*Pi)}; //              // rotation angle array

    for (int gly_id = 0; gly_id< nglyphs; gly_id++ )
    {   
        int glyph_cnt = gly_id * n_coords/3;
        int index = 3*gly_id;

        kvs::Vector3f dir(direction.at(index),
                          direction.at(index+1),
                          direction.at(index+2));
        // translate & scaling option
        const kvs::Vector3f trans_direction = dir.normalized();
        const kvs::Vec3 v = trans_direction.normalized();
        const kvs::Vec3 c = DefaultDirection.cross( v );
        const float d = DefaultDirection.dot( v );
        const float s = static_cast<float>( std::sqrt( ( 1.0 + d ) * 2.0 ) );
        const kvs::Quaternion q( c.x()/s, c.y()/s, c.z()/s, s/2.0f );
        const kvs::Mat3 trans_scale = q.toMatrix();
        const kvs::Vector3f trans_position(position.data()+ index);
        const kvs::Xform xform_trans( trans_position, BaseClass::scale() * m_size, trans_scale );
#ifdef DEBUG        
        std::cout << "dir = " << dir.x() << ", "<< dir.y() << ", "<<dir.z() << std::endl;
        std::cout << "trans_direction = " << trans_direction.x() << ", "<< trans_direction.y() << ", "<< trans_direction.z() << std::endl;
        std::cout << "v = " << v.x() << ", "<< v.y() << ", "<< v.z() << std::endl;
        std::cout << "c = " << c.x() << ", "<< c.y() << ", "<< c.z() << std::endl;
        std::cout << "q = " << q.x() << ", "<< q.y() << ", "<< q.z() << ", " << q.w()<< std::endl;
#endif

        // rotate option
        float x= V.x();
        float y= V.y();
        float z= V.z();
//        float cos=std::cos(R[gly_id]);
//        float sin=std::sin(R[gly_id]);
        float cos=std::cos(R);
        float sin=std::sin(R);

        kvs::Mat4 rot( x*x*(1.0f-cos)+cos, x*y*(1.0f-cos)-z*sin, x*z*(1.0f-cos)+y*sin, 0,
                       y*x*(1.0f-cos)*-z*sin, y*y*(1.0f-cos)+cos,y*z*(1.0f-cos)-x*sin, 0,
                       x*z*(1.0f-cos)*-y*sin, y*z*(1.0f-cos)+x*sin, z*z*(1.0f-cos)+cos,0,
                       0, 0, 0, 1); 
        kvs::Xform xform_rot(rot);
        
#ifdef DEBUG
        for(int i =0; i<4 ;i++)
        {
            std::cout << "xform_trans.toMatrix[0]["<<i<<"]=" <<xform_trans.toMatrix()[0][i] << std::endl;
            std::cout << "xform_trans.toMatrix[1]["<<i<<"]=" <<xform_trans.toMatrix()[1][i] << std::endl;
            std::cout << "xform_trans.toMatrix[2]["<<i<<"]=" <<xform_trans.toMatrix()[2][i] << std::endl;
            std::cout << "xform_trans.toMatrix[3]["<<i<<"]=" <<xform_trans.toMatrix()[3][i] << std::endl;
        }
#endif
        for (int id = 0; id < n_coords; id +=3)
        {
            kvs::Vec3 co(m_coords.at(id), m_coords.at(id+1), m_coords.at(id+2)); 
            kvs::Vec3 co_normal(m_normal.at(id), m_normal.at(id+1), m_normal.at(id+2)); 

            //rotate
            co        = co*xform_rot.rotation();
            co_normal = co_normal*xform_rot.rotation(); 

            //scale & translate
//            co        = co*xform_trans.rotation();
//            co_normal = co_normal*xform_trans.rotation(); 
            co        = xform_trans.transform(co);
            co_normal = xform_trans.transformNormal(co_normal);

            // set global value
            gl_coords.at(gly_id*n_coords+id  )= co[0]; 
            gl_coords.at(gly_id*n_coords+id+1)= co[1]; 
            gl_coords.at(gly_id*n_coords+id+2)= co[2]; 
            gl_normal.at(gly_id*n_coords+id  )= co_normal[0]; 
            gl_normal.at(gly_id*n_coords+id+1)= co_normal[1]; 
            gl_normal.at(gly_id*n_coords+id+2)= co_normal[2]; 
        }
        
        for (int id = 0;id < n_connection; id+=3 )
        {
            gl_connection.at(gly_id*n_connection+id  )= m_connection[id+0]+glyph_cnt; 
            gl_connection.at(gly_id*n_connection+id+1)= m_connection[id+1]+glyph_cnt; 
            gl_connection.at(gly_id*n_connection+id+2)= m_connection[id+2]+glyph_cnt; 
        }
    }

    m_coords.allocate(gl_coords.size());
    m_normal.allocate(gl_normal.size());
    m_connection.allocate(gl_connection.size()); 

    m_coords = gl_coords;
    m_normal = gl_normal;
    m_connection = gl_connection;

#ifdef DEBUG
    std::cout << "object->coord()  = " << float(gl_coords[0]) << ", " <<  float(gl_coords[1]) << ", " << float(gl_coords[2]) << ", "  << 
        float(gl_coords[402]) << ", " <<  float(gl_coords[403]) << ", " << float(gl_coords[404]) << ", " <<
        float(gl_coords[564]) << ", " <<  float(gl_coords[565]) << ", " << float(gl_coords[566]) << ", " << std::endl;
    std::cout << "connection[id]  = " << float(gl_connection[882]) << ", " <<  float(gl_connection[883]) << ", " << float(gl_connection[884]) << ", "  << 
        float(gl_connection[1761]) << ", " <<  float(gl_connection[1762]) << ", " << float(gl_connection[1763]) << ", " << std::endl;
    std::cout  <<  __PRETTY_FUNCTION__ <<" : "<<__LINE__ << std::endl;
#endif

    // set polygon
    SuperClass::setCoords( m_coords  ); 
    SuperClass::setConnections( m_connection );
    SuperClass::setColor( m_color );
    SuperClass::setNormals( m_normal );
    SuperClass::setOpacity( int(m_opacity) );
    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::PolygonColor );
}

/*===========================================================================*/
/**
 *  @brief  setting arrow glyph polygon.
 */
/*===========================================================================*/
//void ArrowGlyph::setPolygon()
//{
////    変換後に格納
//    SuperClass::setCoords( m_coords  ); 
//    SuperClass::setConnections( m_connection );
//    SuperClass::setColor( m_color );
//    SuperClass::setNormals( m_normal );
//    SuperClass::setOpacity( int(m_opacity) );
//    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
//    SuperClass::setColorType( kvs::PolygonObject::PolygonColor );
//    std::cout  <<  __PRETTY_FUNCTION__ <<" : "<<__LINE__ << std::endl;
////#endif
//
//}

/*===========================================================================*/
/**
 *  @brief  Initialize OpenGL properties for rendering arrow glyph.
 */
/*===========================================================================*/
void ArrowGlyph::initialize()
{
    kvs::OpenGL::SetBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    kvs::OpenGL::SetPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    kvs::OpenGL::SetShadeModel( GL_SMOOTH );
    kvs::OpenGL::SetColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
    kvs::OpenGL::SetLightModel( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );

    kvs::OpenGL::Disable( GL_LINE_SMOOTH );
    kvs::OpenGL::Enable( GL_BLEND );
    kvs::OpenGL::Enable( GL_COLOR_MATERIAL );
    if ( m_type == ArrowGlyph::LineArrow )
    {
        kvs::OpenGL::Disable( GL_NORMALIZE );
        kvs::OpenGL::Disable( GL_LIGHTING );
    }
    else
    {
        if ( !BaseClass::isEnabledShading() )
        {
            kvs::OpenGL::Disable( GL_NORMALIZE );
            kvs::OpenGL::Disable( GL_LIGHTING );
        }
        else
        {
            kvs::OpenGL::Enable( GL_NORMALIZE );
            kvs::OpenGL::Enable( GL_LIGHTING );
        }
    }
}

} // end of namespace kvs
