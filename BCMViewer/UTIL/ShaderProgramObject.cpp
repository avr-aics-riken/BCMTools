/*
 *  BCMViewer
 *
 *  Copyright 2012 SGI Japan, Ltd.
 *
 */

#include "ShaderProgramObject.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glext.h>
#endif

using namespace std;

namespace
{

// シェーダの情報を表示する
void printShaderInfoLog(GLuint shader)
{
	GLsizei bufSize;
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH , &bufSize);
	if (bufSize > 1)
	{
		GLchar *infoLog;
		infoLog = (GLchar *)malloc(bufSize);
		if (infoLog != NULL)
		{
			GLsizei length;
			glGetShaderInfoLog(shader, bufSize, &length, infoLog);
			cout << "InfoLog:" << endl << infoLog << endl << endl;
			free(infoLog);
		}
		else
		{
			cout << "Could not allocate InfoLog buffer." << endl;
		}
	}
}

// プログラムの情報を表示する
void printProgramInfoLog(GLuint program)
{
	GLsizei bufSize;
	glGetProgramiv(program, GL_INFO_LOG_LENGTH , &bufSize);
	if (bufSize > 1)
	{
		GLchar *infoLog;
		infoLog = (GLchar *)malloc(bufSize);
		if (infoLog != NULL)
		{
			GLsizei length;
			glGetProgramInfoLog(program, bufSize, &length, infoLog);
			cout << "InfoLog:" << endl << infoLog << endl << endl;
			free(infoLog);
		}
		else
		{
			cout << "Could not allocate InfoLog buffer." << endl;
		}
	}
}

};// namespace


ShaderObject::ShaderObject()
{
	m_shader = 0;
}

ShaderObject::~ShaderObject()
{
	Release();
}


bool ShaderObject::LoadFromFile(const std::string& filename, SHADERTYPE shaderType)
{
	std::ifstream fin(filename.c_str(), ios::binary);
	if (fin.fail())
	{
		cout << "cannot open shader file: " << filename << endl;
		return false;
	}
	std::ostringstream str_out;
	str_out << fin.rdbuf();
	std::string fileBuf = str_out.str();
	fin.close();
	return LoadFromMemory(fileBuf, shaderType);
}

bool ShaderObject::LoadFromMemory(const std::string& programSource, SHADERTYPE shaderType)
{
	const std::string prgSource(programSource);
	m_source = prgSource;

	const char* s = prgSource.c_str();
	int l = static_cast<int>(prgSource.length());
	
	m_shader = glCreateShader(shaderType);
	glShaderSource( m_shader, 1, &s, &l );
	if ( glGetError() != GL_NO_ERROR )
	{
		cout << "cannot set shader source: " << prgSource << endl;
		return false;
	}

	// compile
	int compiled = 0;
	glCompileShader(m_shader);
	glGetShaderiv(m_shader, GL_COMPILE_STATUS, &compiled);
	printShaderInfoLog(m_shader);
	if (!compiled)
	{
		cout << "Compile is failed" << endl;
		return false;
	}

	return true;
}

void ShaderObject::Release()
{
	if (m_shader)
	{
		glDeleteShader(m_shader);
		m_shader = 0;
	}
}


ProgramObject::ProgramObject()
{
	m_oldProgram = 0;
	m_program = 0;
	m_binding = false;
}

ProgramObject::ProgramObject(const ShaderObject& vertexShader, const ShaderObject& fragmentShader)
{
	m_oldProgram = 0;
	m_program = 0;
	m_binding = false;
	Link(vertexShader, fragmentShader);
}

bool ProgramObject::Link(const ShaderObject& vertexShader, const ShaderObject& fragmentShader)
{
	unsigned int glProgram = glCreateProgram();
	glAttachShader(glProgram, vertexShader.GetShader());
	glAttachShader(glProgram, fragmentShader.GetShader());

	/* シェーダプログラムのリンク */
	glLinkProgram(glProgram);
	GLint linked;
	glGetProgramiv(glProgram, GL_LINK_STATUS, &linked);
	printProgramInfoLog(glProgram);
	if (linked == GL_FALSE)
	{
		cout << "Link error." << endl;
		return false;
	}
	m_program = glProgram;
	return true;
}

void ProgramObject::Bind()
{
	glGetIntegerv(GL_CURRENT_PROGRAM, &m_oldProgram); 	
	glUseProgram(m_program);
	m_binding = true;
}

void ProgramObject::Unbind()
{
	glUseProgram(m_oldProgram);
	m_binding = false;
}

void ProgramObject::Release()
{
	if (m_program)
	{
		glDeleteProgram(m_program);
		m_program = 0;
	}
}

// 1i - 4i
void ProgramObject::SetUniform(const char* name, const int i0)
{
	if (m_binding)
		glUniform1i(glGetUniformLocation(m_program, name), i0);
}

void ProgramObject::SetUniform(const char* name, const int i0, const int i1)
{
	if (m_binding)
		glUniform2i(glGetUniformLocation(m_program, name), i0, i1);
}

void ProgramObject::SetUniform(const char* name, const int i0, const int i1, const int i2)
{
	if (m_binding)
		glUniform3i(glGetUniformLocation(m_program, name), i0, i1, i2);
}

void ProgramObject::SetUniform(const char* name, const int i0, const int i1, const int i2, const int i3)
{
	if (m_binding)
		glUniform4i(glGetUniformLocation(m_program, name), i0, i1, i2, i3);
}
/*
void ProgramObject::SetUniform(const char* name, const int num, const int* i_array)
{
	void (APIENTRYP glUniFuncs[]) (GLint, GLsizei, const GLint*) = 
	{
		0,
		glUniform1iv,
		glUniform2iv,
		glUniform3iv,
		glUniform4iv
	};
	if (num <= 0 || num >= 5)
		return;

	if (m_binding)
		glUniFuncs[num](glGetUniformLocation(m_program, name), num, i_array);
}
*/

// 1f - 4f
void ProgramObject::SetUniform(const char* name, const float f0)
{
	if (m_binding)
		glUniform1f(glGetUniformLocation(m_program, name), f0);
}

void ProgramObject::SetUniform(const char* name, const float f0, const float f1)
{
	if (m_binding)
		glUniform2f(glGetUniformLocation(m_program, name), f0, f1);
}

void ProgramObject::SetUniform(const char* name, const float f0, const float f1, const float f2)
{
	if (m_binding)
		glUniform3f(glGetUniformLocation(m_program, name), f0, f1, f2);
}

void ProgramObject::SetUniform(const char* name, const float f0, const float f1, const float f2, const float f3)
{
	if (m_binding)
		glUniform4f(glGetUniformLocation(m_program, name), f0, f1, f2, f3);
}

/*
void ProgramObject::SetUniform(const char* name, const int num, const float* f_array)
{
	void (APIENTRYP glUniFuncs[]) (GLint, GLsizei, const GLfloat*) = 
	{
		0,
		glUniform1fv,
		glUniform2fv,
		glUniform3fv,
		glUniform4fv
	};
	if (num <= 0 || num >= 5)
		return;

	if (m_binding)
		glUniFuncs[num](glGetUniformLocation(m_program, name), num, f_array);
}
*/

// matrix
void ProgramObject::SetUniformMatrix2x2(const char* name, const int count, const bool transpose, const float* val)
{
	if (m_binding)
		glUniformMatrix2fv(glGetUniformLocation(m_program, name), count, transpose, val);
}
void ProgramObject::SetUniformMatrix3x3(const char* name, const int count, const bool transpose, const float* val)
{
	if (m_binding)
		glUniformMatrix3fv(glGetUniformLocation(m_program, name), count, transpose, val);
}
void ProgramObject::SetUniformMatrix4x4(const char* name, const int count, const bool transpose, const float* val)
{
	if (m_binding)
		glUniformMatrix4fv(glGetUniformLocation(m_program, name), count, transpose, val);
}


