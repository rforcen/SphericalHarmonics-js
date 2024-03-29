// mesh object argument { coords:, normals:, colors:, textures: }
// of float32array x 3 and textures x 2

function main(mesh) {
  const canvas = document.querySelector('#canvas')
  const gl = canvas.getContext('webgl')

  // If we don't have a GL context, give up now
  if (!gl) {
    alert('Unable to initialize WebGL. Your browser or machine may not support it.')
    return
  }

  // Vertex shader program
  const vsSource = `
    attribute vec4 aVertexPosition, aVertexColor, aVertexNormal; // vertex, color, normal
    uniform mat4 uModelViewMatrix, uProjectionMatrix; // view mtx

    varying lowp vec4 vColor, vNormal, fragPos; // -> fragment

    void main(void) {
      vColor = aVertexColor; 
      vNormal = aVertexNormal;
      fragPos = uModelViewMatrix * aVertexPosition; 

      gl_Position = uProjectionMatrix * fragPos;
    }
  `

  // Fragment shader program
  const fsSource = `
    varying lowp vec4 vColor, vNormal, fragPos;

    lowp vec4 calc_light() { // vColor, vNormal
      lowp float l=0.8;
      lowp vec4 lightColor = vec4(l, l, l, 1);
      lowp vec4 lightPos = vec4(1, 0.5, 1, 1);

      // ambient
      lowp vec4 ambient = 0.2 * lightColor;

      // Diffuse 
      lowp vec4 lightDir = normalize(lightPos - fragPos);
      lowp float diff = max(dot(vNormal, lightDir), 0.0);
      lowp vec4 diffuse = diff * lightColor;

      return (ambient + diffuse) * vColor;
    }

    void main(void) {
      gl_FragColor = calc_light();
    }
  `

  // Initialize a shader program this is where all the lighting
  // for the vertices and so forth is established.
  const shaderProgram = initShaderProgram(gl, vsSource, fsSource)

  // Collect all the info needed to use the shader program.
  // Look up which attributes our shader program is using
  // for aVertexPosition, aVevrtexColor and also
  // look up uniform locations.
  const programInfo = {
    program: shaderProgram,
    attribLocations: {
      vertexPosition: gl.getAttribLocation(shaderProgram, 'aVertexPosition'),
      vertexColor: gl.getAttribLocation(shaderProgram, 'aVertexColor'),
      vertexNormal: gl.getAttribLocation(shaderProgram, 'aVertexNormal'),
    },
    uniformLocations: {
      projectionMatrix: gl.getUniformLocation(shaderProgram, 'uProjectionMatrix'),
      modelViewMatrix: gl.getUniformLocation(shaderProgram, 'uModelViewMatrix'),
    },
  }

  // Here's where we call the routine that builds all the
  // objects we'll be drawing.
  const buffers = initBuffers(gl, mesh)

  var then = 0

  // Draw the scene repeatedly
  function render(now) {
    now *= delay // convert to seconds
    const deltaTime = now - then
    then = now

    drawScene(gl, programInfo, buffers, deltaTime)

    requestAnimationFrame(render)
  }


  requestAnimationFrame(render)
}

//
// initBuffers
//

function initBuffers(gl, mesh) {

  function create_buffer(v) {
    const buff = gl.createBuffer()
    gl.bindBuffer(gl.ARRAY_BUFFER, buff)
    gl.bufferData(gl.ARRAY_BUFFER, v, gl.STATIC_DRAW)
    return buff
  }

  return {
    position: create_buffer(mesh.coords),
    normal: create_buffer(mesh.normals),
    color: create_buffer(mesh.colors)
  }
}

//
// Draw the scene.
//
function drawScene(gl, programInfo, buffers, deltaTime) {
  function init() {
    gl.clearColor(0.0, 0.0, 0.0, 1.0) // Clear to black, fully opaque
    gl.clearDepth(1.0) // Clear everything
    gl.enable(gl.DEPTH_TEST) // Enable depth testing
    gl.depthFunc(gl.LEQUAL) // Near things obscure far things
  }

  init()

  // Clear the canvas before we start drawing on it.
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

  // Create a perspective matrix, a special matrix that is
  // used to simulate the distortion of perspective in a camera.
  // Our field of view is 45 degrees, with a width/height
  // ratio that matches the display size of the canvas
  // and we only want to see objects between 0.1 units
  // and 100 units away from the camera.

  const fieldOfView = 45 * Math.PI / 180 // in radians
  const aspect = gl.canvas.clientWidth / gl.canvas.clientHeight
  const zNear = 0.1
  const zFar = 100.0
  const projectionMatrix = mat4.create()

  // note: glmatrix.js always has the first argument
  // as the destination to receive the result.
  mat4.perspective(projectionMatrix, fieldOfView, aspect, zNear, zFar)

  // Set the drawing position to the "identity" point, which is
  // the center of the scene.
  const modelViewMatrix = mat4.create()

  // Now move the drawing position a bit to where we want to
  // start drawing the square.

  const zoom = -3.5
  mat4.translate(modelViewMatrix, modelViewMatrix, [0, 0, zoom]) // amount to translate
  mat4.rotate(modelViewMatrix, modelViewMatrix, squareRotation, [0.5, 0.4, 0.1]) // axis to rotate around


  function setAttrib(btype, attrib, numComponents) { // setAttrib(buffers.position, programInfo.attribLocations.vertexPosition, 3)
    const type = gl.FLOAT
    const normalize = false
    const stride = 0
    const offset = 0
    gl.bindBuffer(gl.ARRAY_BUFFER, btype)
    gl.vertexAttribPointer(
      attrib,
      numComponents,
      type,
      normalize,
      stride,
      offset)
    gl.enableVertexAttribArray(attrib)
  }

  setAttrib(buffers.position, programInfo.attribLocations.vertexPosition, 3)
  setAttrib(buffers.color, programInfo.attribLocations.vertexColor, 3)
  setAttrib(buffers.normal, programInfo.attribLocations.vertexNormal, 3)

  // Tell WebGL to use our program when drawing
  gl.useProgram(programInfo.program)

  // Set the shader uniforms
  gl.uniformMatrix4fv(programInfo.uniformLocations.projectionMatrix, false, projectionMatrix)
  gl.uniformMatrix4fv(programInfo.uniformLocations.modelViewMatrix, false, modelViewMatrix)

  gl.drawArrays(gl.TRIANGLE_FAN, 0, sh.n_vertex)

  // Update the rotation for the next draw
  squareRotation += deltaTime
}

//
// Initialize a shader program, so WebGL knows how to draw our data
//
function initShaderProgram(gl, vsSource, fsSource) {
  const vertexShader = loadShader(gl, gl.VERTEX_SHADER, vsSource)
  const fragmentShader = loadShader(gl, gl.FRAGMENT_SHADER, fsSource)

  // Create the shader program

  const shaderProgram = gl.createProgram()
  gl.attachShader(shaderProgram, vertexShader)
  gl.attachShader(shaderProgram, fragmentShader)
  gl.linkProgram(shaderProgram)

  // If creating the shader program failed, alert

  if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
    alert('Unable to initialize the shader program: ' + gl.getProgramInfoLog(shaderProgram))
    return null
  }

  return shaderProgram
}

//
// creates a shader of the given type, uploads the source and
// compiles it.
//
function loadShader(gl, type, source) {
  const shader = gl.createShader(type)

  // Send the source to the shader object
  gl.shaderSource(shader, source)

  // Compile the shader program

  gl.compileShader(shader)

  // See if it compiled successfully

  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    alert('An error occurred compiling the shaders: ' + gl.getShaderInfoLog(shader))
    gl.deleteShader(shader)
    return null
  }

  return shader
}