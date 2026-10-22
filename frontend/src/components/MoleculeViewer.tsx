import { useEffect, useMemo, useRef, useState } from 'react'
import { Expand, LoaderCircle, Rotate3D, X } from 'lucide-react'
import * as THREE from 'three'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'
import { getStructure } from '../lib/api'
import type { StructureData } from '../lib/types'

const elementStyle: Record<string, { color: number; radius: number }> = {
  H: { color: 0xe8f1ed, radius: 0.24 },
  C: { color: 0x394841, radius: 0.36 },
  O: { color: 0xf06f65, radius: 0.34 },
  Na: { color: 0x8c82f5, radius: 0.48 },
}

interface Props {
  jobId: string | null
  structureReference: string | null
  label: string
}

export function MoleculeViewer({ jobId, structureReference, label }: Props) {
  const mountRef = useRef<HTMLDivElement>(null)
  const [structure, setStructure] = useState<StructureData | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const [expanded, setExpanded] = useState(false)

  useEffect(() => {
    setExpanded(false)
    setStructure(null)
    setError(null)
    if (!jobId || !structureReference) return
    let cancelled = false
    setLoading(true)
    getStructure(jobId, structureReference)
      .then((value) => { if (!cancelled) setStructure(value) })
      .catch((caught) => { if (!cancelled) setError(caught instanceof Error ? caught.message : 'Could not load geometry') })
      .finally(() => { if (!cancelled) setLoading(false) })
    return () => { cancelled = true }
  }, [jobId, structureReference])

  useEffect(() => {
    const mount = mountRef.current
    if (!mount || !structure) return
    mount.replaceChildren()

    const scene = new THREE.Scene()
    scene.fog = new THREE.FogExp2(0x07100d, 0.018)
    const camera = new THREE.PerspectiveCamera(38, 1, 0.1, 200)
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true })
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2))
    renderer.outputColorSpace = THREE.SRGBColorSpace
    renderer.setClearColor(0x07100d, 1)
    mount.appendChild(renderer.domElement)

    const group = new THREE.Group()
    scene.add(group)
    scene.add(new THREE.HemisphereLight(0xd7fff0, 0x132019, 2.2))
    const key = new THREE.DirectionalLight(0xffffff, 2.8)
    key.position.set(5, 8, 6)
    scene.add(key)
    const rim = new THREE.DirectionalLight(0x73efb7, 2.0)
    rim.position.set(-6, -2, 5)
    scene.add(rim)

    const atomPositions = new Map<number, THREE.Vector3>()
    for (const atom of structure.atoms) {
      const position = new THREE.Vector3(...atom.position)
      atomPositions.set(atom.index, position)
      const style = elementStyle[atom.element] ?? { color: 0x74b5a0, radius: 0.34 }
      const geometry = new THREE.SphereGeometry(style.radius, 28, 20)
      const material = new THREE.MeshStandardMaterial({
        color: style.color,
        roughness: 0.28,
        metalness: atom.element === 'Na' ? 0.5 : 0.08,
      })
      const mesh = new THREE.Mesh(geometry, material)
      mesh.position.copy(position)
      mesh.userData = { element: atom.element, index: atom.index }
      group.add(mesh)
    }

    for (const bond of structure.bonds) {
      const start = atomPositions.get(bond.source)
      const end = atomPositions.get(bond.target)
      if (!start || !end) continue
      const direction = new THREE.Vector3().subVectors(end, start)
      const length = direction.length()
      const geometry = new THREE.CylinderGeometry(0.075, 0.075, length, 12)
      const material = new THREE.MeshStandardMaterial({ color: 0x8da69b, roughness: 0.42 })
      const cylinder = new THREE.Mesh(geometry, material)
      cylinder.position.copy(start).add(end).multiplyScalar(0.5)
      cylinder.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize())
      group.add(cylinder)
    }

    const bounds = new THREE.Box3().setFromObject(group)
    const center = bounds.getCenter(new THREE.Vector3())
    const size = Math.max(bounds.getSize(new THREE.Vector3()).length(), 2.5)
    group.position.sub(center)
    camera.position.set(size * 0.72, size * 0.52, size * 1.05)
    camera.lookAt(0, 0, 0)

    const controls = new OrbitControls(camera, renderer.domElement)
    controls.enableDamping = true
    controls.dampingFactor = 0.06
    controls.autoRotate = true
    controls.autoRotateSpeed = 0.65
    controls.minDistance = 2
    controls.maxDistance = Math.max(size * 4, 12)

    const resize = () => {
      const width = Math.max(mount.clientWidth, 1)
      const height = Math.max(mount.clientHeight, 1)
      renderer.setSize(width, height, false)
      camera.aspect = width / height
      camera.updateProjectionMatrix()
    }
    const observer = new ResizeObserver(resize)
    observer.observe(mount)
    resize()

    let frame = 0
    const animate = () => {
      frame = requestAnimationFrame(animate)
      controls.update()
      renderer.render(scene, camera)
    }
    animate()

    return () => {
      cancelAnimationFrame(frame)
      observer.disconnect()
      controls.dispose()
      group.traverse((child) => {
        if (child instanceof THREE.Mesh) {
          child.geometry.dispose()
          const materials = Array.isArray(child.material) ? child.material : [child.material]
          materials.forEach((material) => material.dispose())
        }
      })
      renderer.dispose()
      renderer.domElement.remove()
    }
  }, [structure, expanded])

  const elements = useMemo(
    () => [...new Set(structure?.atoms.map((atom) => atom.element) ?? [])],
    [structure],
  )

  return (
    <section className={`molecule-viewer ${expanded ? 'expanded' : ''}`}>
      {expanded && <button className="viewer-backdrop" aria-label="Close expanded viewer" onClick={() => setExpanded(false)} />}
      <div className="viewer-card">
        <div className="viewer-heading">
          <div><span>REPRESENTATIVE GEOMETRY</span><strong>{label}</strong></div>
          <div className="viewer-actions">
            <span><Rotate3D size={12} /> drag · zoom</span>
            <button onClick={() => setExpanded(!expanded)} aria-label={expanded ? 'Close expanded molecule viewer' : 'Expand molecule viewer'}>{expanded ? <X size={15} /> : <Expand size={15} />}</button>
          </div>
        </div>
        <div className="molecule-canvas" ref={mountRef}>
          {loading && <div className="viewer-state"><LoaderCircle size={19} className="spin" /> Loading saved geometry…</div>}
          {error && <div className="viewer-state error">{error}</div>}
          {!loading && !error && !structure && <div className="viewer-state">No representative structure saved</div>}
        </div>
        {structure && <div className="element-legend">
          {elements.map((element) => <span key={element}><i style={{ background: `#${(elementStyle[element]?.color ?? 0x74b5a0).toString(16).padStart(6, '0')}` }} />{element}</span>)}
          <em>{structure.atoms.length} atoms · {structure.bonds.length} detected bonds</em>
        </div>}
        {expanded && structure && <p className="geometry-note">{structure.note}</p>}
      </div>
    </section>
  )
}
