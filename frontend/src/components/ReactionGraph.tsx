import { useEffect, useMemo } from 'react'
import {
  Background,
  BaseEdge,
  EdgeLabelRenderer,
  Controls,
  Handle,
  MarkerType,
  Position,
  ReactFlow,
  useNodesState,
  type Edge,
  type EdgeProps,
  type Node,
  type NodeProps,
} from '@xyflow/react'
import type { GraphEdge, GraphNode, Selection } from '../lib/types'

function ChemistryNode({ data, selected }: NodeProps<Node<{ record: GraphNode }>>) {
  const record = data.record
  return (
    <div className={`chem-node ${record.role} ${record.interpretation_quality === 'uncertain' ? 'uncertain' : ''} ${selected ? 'selected' : ''}`}>
      <Handle id="in" type="target" position={Position.Left} />
      <Handle id="return-in" type="target" position={Position.Top} />
      <span className="node-role">{record.role === 'endpoint' ? 'search endpoint' : record.role}</span>
      <strong>{record.label}</strong>
      {record.interpretation_quality === 'uncertain' && <span className="node-warning">unresolved connectivity</span>}
      <div className="node-metrics">
        <span>ΔE {record.relative_energy_ev.toFixed(2)} eV</span>
        <span>seen {record.runs_seen} · stops {record.terminal_runs}</span>
      </div>
      <Handle id="out" type="source" position={Position.Right} />
      <Handle id="return-out" type="source" position={Position.Top} style={{ left: '65%' }} />
    </div>
  )
}

const nodeTypes = { chemistry: ChemistryNode }

function ReturnEdge({ id, sourceX, sourceY, targetX, targetY, markerEnd, style, label }: EdgeProps) {
  const top = Math.min(sourceY, targetY) - 170
  const path = `M ${sourceX} ${sourceY} C ${sourceX} ${top}, ${targetX} ${top}, ${targetX} ${targetY}`
  return <><BaseEdge id={id} path={path} markerEnd={markerEnd} style={style} />
    <EdgeLabelRenderer><div className="return-edge-label nodrag nopan" style={{ transform: `translate(-50%, -50%) translate(${(sourceX + targetX) / 2}px, ${(sourceY + targetY) / 8 + top * .75}px)` }}>{label}</div></EdgeLabelRenderer>
  </>
}
const edgeTypes = { returnPath: ReturnEdge }

function calculateDepths(records: GraphNode[], edges: GraphEdge[]): Map<string, number> {
  const depth = new Map<string, number>()
  const roots = records.filter((record) => record.starting_runs > 0 || record.role === 'reactant')
  const queue = (roots.length ? roots : records.slice(0, 1)).map((record) => record.id)
  queue.forEach((id) => depth.set(id, 0))
  const outgoing = new Map<string, string[]>()
  edges.forEach((edge) => outgoing.set(edge.source, [...(outgoing.get(edge.source) ?? []), edge.target]))
  while (queue.length) {
    const source = queue.shift()!
    const nextDepth = (depth.get(source) ?? 0) + 1
    for (const target of outgoing.get(source) ?? []) {
      if (depth.has(target)) continue
      depth.set(target, nextDepth)
      queue.push(target)
    }
  }
  const fallbackDepth = Math.max(0, ...depth.values()) + 1
  records.forEach((record) => { if (!depth.has(record.id)) depth.set(record.id, fallbackDepth) })
  return depth
}

interface Props {
  nodes: GraphNode[]
  edges: GraphEdge[]
  minimumFrequency: number
  totalRuns: number
  onSelect: (selection: Selection) => void
}

export function ReactionGraph({ nodes: records, edges: recordsEdges, minimumFrequency, totalRuns, onSelect }: Props) {
  const visibleEdges = useMemo(
    () => recordsEdges.filter((edge) => edge.run_frequency * 100 >= minimumFrequency),
    [recordsEdges, minimumFrequency],
  )
  const visibleIds = useMemo(() => new Set(visibleEdges.flatMap((edge) => [edge.source, edge.target])), [visibleEdges])
  const visibleRecords = useMemo(
    () => minimumFrequency === 0 ? records : records.filter((node) => visibleIds.has(node.id) || node.starting_runs > 0 || node.role === 'reactant'),
    [records, visibleIds, minimumFrequency],
  )
  const computedNodes = useMemo(() => {
    const depths = calculateDepths(visibleRecords, visibleEdges)
    const rows = new Map<number, number>()
    return visibleRecords.map((record): Node<{ record: GraphNode }> => {
      const column = depths.get(record.id) ?? 0
      const row = rows.get(column) ?? 0
      rows.set(column, row + 1)
      return {
        id: record.id,
        type: 'chemistry',
        position: { x: column * 420 + 35, y: row * 240 + (column % 2) * 145 + 180 },
        data: { record },
      }
    })
  }, [visibleRecords, visibleEdges])
  const [flowNodes, setFlowNodes, onNodesChange] = useNodesState(computedNodes)
  useEffect(() => {
    setFlowNodes((current) => computedNodes.map((node) => {
      const existing = current.find((candidate) => candidate.id === node.id)
      return existing ? { ...node, position: existing.position } : node
    }))
  }, [computedNodes, setFlowNodes])

  const flowEdges: Edge[] = visibleEdges.map((edge) => ({
    id: edge.id,
    source: edge.source,
    target: edge.target,
    sourceHandle: edge.source > edge.target ? 'return-out' : 'out',
    targetHandle: edge.source > edge.target ? 'return-in' : 'in',
    type: edge.source > edge.target ? 'returnPath' : 'smoothstep',
    animated: edge.delta_energy_ev < 0,
    label: `seen ${Math.round(edge.run_frequency * 100)}% (${edge.run_count}/${totalRuns}) · ΔE ${edge.delta_energy_ev >= 0 ? '+' : ''}${edge.delta_energy_ev.toFixed(2)} eV`,
    labelStyle: { fill: '#d2e0da', fontSize: 12, fontFamily: 'ui-monospace, monospace' },
    labelBgStyle: { fill: '#0b1612', fillOpacity: 0.92 },
    style: { stroke: '#80e3b5', strokeWidth: 1.5 + edge.run_frequency * 3 },
    markerEnd: { type: MarkerType.ArrowClosed, color: '#80e3b5' },
    data: { record: edge },
  }))

  return (
    <ReactFlow
      nodes={flowNodes}
      edges={flowEdges}
      onNodesChange={onNodesChange}
      nodesDraggable
      nodesConnectable={false}
      nodeTypes={nodeTypes}
      edgeTypes={edgeTypes}
      fitView
      fitViewOptions={{ padding: 0.18 }}
      minZoom={0.35}
      maxZoom={1.7}
      onNodeClick={(_, node) => onSelect({ kind: 'node', value: node.data.record })}
      onEdgeClick={(_, edge) => onSelect({ kind: 'edge', value: (edge.data as { record: GraphEdge }).record })}
      proOptions={{ hideAttribution: true }}
    >
      <Background color="#1c342b" gap={22} size={1} />
      <Controls showInteractive={false} />
    </ReactFlow>
  )
}
