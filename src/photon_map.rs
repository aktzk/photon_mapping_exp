use glm::DVec3;
use ordered_float::OrderedFloat;
use std::collections::BinaryHeap;
use std::sync::Arc;

pub struct Photon {
    pub pos: DVec3,
    pub power: DVec3,
    pub incoming_dir: DVec3,
}

impl Photon {
    pub fn new(pos: DVec3, power: DVec3, incoming_dir: DVec3) -> Self {
        Self {
            pos,
            power,
            incoming_dir,
        }
    }
}

pub struct PhotonMap {
    root_node: Option<Box<TreeNode>>,
}

impl PhotonMap {
    pub fn new(photons: &mut Vec<Arc<Photon>>) -> Self {
        let end = photons.len();
        Self {
            root_node: construct_kd_tree(photons, 0, end, 0),
        }
    }

    pub fn collect_photons(&self, arg: &PhotonCollectArg) -> NearestPhotons {
        let mut result = NearestPhotons::new();
        if let Some(node) = &self.root_node {
            self.collect_photons_inner(
                &mut result,
                &node,
                arg.max_distance2,
                arg.max_collect_num,
                &arg.center_pos,
                &arg.normal,
            );
        }
        result
    }

    fn collect_photons_inner(
        &self,
        result: &mut NearestPhotons,
        node: &TreeNode,
        max_distance2: f64,
        max_collect_num: usize,
        center_pos: &glm::DVec3,
        normal: &glm::DVec3,
    ) {
        let photon = &node.photon;

        let delta = match node.axis {
            NodeAxis::X => center_pos.x - photon.pos.x,
            NodeAxis::Y => center_pos.y - photon.pos.y,
            NodeAxis::Z => center_pos.z - photon.pos.z,
        };

        let dir = photon.pos - center_pos;
        let distance2 = glm::length2(&dir);
        let normal_dot_dir = normal.dot(&(dir / distance2.sqrt()));

        // 楕円状の近傍を集める
        if distance2 < max_distance2 && normal_dot_dir.abs() < max_distance2 * 0.01 {
            result
                .queue
                .push(NeighborPhoton::new(photon.clone(), distance2));
            if result.queue.len() > max_collect_num {
                let _ = result.queue.pop();
            }
        }
        let new_max_distance2 = result
            .queue
            .peek()
            .and_then(|p| Some(p.distance2.into_inner()))
            .unwrap_or(max_distance2);
        if delta > 0.0 {
            if let Some(right) = node.right_node.as_ref() {
                self.collect_photons_inner(
                    result,
                    right,
                    new_max_distance2,
                    max_collect_num,
                    center_pos,
                    normal,
                );
            }
            if delta * delta < new_max_distance2 {
                if let Some(left) = node.left_node.as_ref() {
                    self.collect_photons_inner(
                        result,
                        left,
                        new_max_distance2,
                        max_collect_num,
                        center_pos,
                        normal,
                    );
                }
            }
        } else {
            if let Some(left) = node.left_node.as_ref() {
                self.collect_photons_inner(
                    result,
                    left,
                    new_max_distance2,
                    max_collect_num,
                    center_pos,
                    normal,
                );
            }
            if delta * delta < new_max_distance2 {
                if let Some(right) = node.right_node.as_ref() {
                    self.collect_photons_inner(
                        result,
                        right,
                        new_max_distance2,
                        max_collect_num,
                        center_pos,
                        normal,
                    );
                }
            }
        }
    }
}

enum NodeAxis {
    X,
    Y,
    Z,
}
struct TreeNode {
    pub photon: Arc<Photon>,
    pub axis: NodeAxis,
    pub left_node: Option<Box<TreeNode>>,
    pub right_node: Option<Box<TreeNode>>,
}

impl TreeNode {
    pub fn new(photon: Arc<Photon>, axis: NodeAxis) -> Self {
        Self {
            photon,
            axis,
            left_node: None,
            right_node: None,
        }
    }
}

// 近傍のフォトンの情報
pub struct NeighborPhoton {
    pub photon: Arc<Photon>,
    pub distance2: OrderedFloat<f64>,
}

impl NeighborPhoton {
    pub fn new(photon: Arc<Photon>, distance2: f64) -> Self {
        Self {
            photon,
            distance2: OrderedFloat(distance2),
        }
    }
}

impl std::cmp::PartialEq for NeighborPhoton {
    fn eq(&self, other: &NeighborPhoton) -> bool {
        self.distance2.eq(&other.distance2)
    }
}
impl std::cmp::Eq for NeighborPhoton {}

impl std::cmp::PartialOrd for NeighborPhoton {
    fn partial_cmp(&self, other: &NeighborPhoton) -> Option<std::cmp::Ordering> {
        self.distance2.partial_cmp(&other.distance2)
    }
}

impl std::cmp::Ord for NeighborPhoton {
    fn cmp(&self, other: &NeighborPhoton) -> std::cmp::Ordering {
        self.distance2.cmp(&other.distance2)
    }
}

// 近傍のフォトン集める
pub struct NearestPhotons {
    pub queue: BinaryHeap<NeighborPhoton>,
}

impl NearestPhotons {
    pub fn new() -> Self {
        Self {
            queue: BinaryHeap::new(),
        }
    }
}

pub struct PhotonCollectArg {
    pub max_distance2: f64,
    pub max_collect_num: usize,
    pub center_pos: glm::DVec3,
    pub normal: glm::DVec3,
}

impl PhotonCollectArg {
    pub fn new(
        max_distance2: f64,
        max_collect_num: usize,
        center_pos: glm::DVec3,
        normal: glm::DVec3,
    ) -> Self {
        Self {
            max_distance2,
            max_collect_num,
            center_pos,
            normal,
        }
    }
}

fn construct_kd_tree(
    photons: &mut Vec<Arc<Photon>>,
    begin: usize,
    end: usize,
    depth: i64,
) -> Option<Box<TreeNode>> {
    if end as i64 - begin as i64 <= 0 {
        return None;
    }
    let axis = match depth % 3 {
        0 => NodeAxis::X,
        1 => NodeAxis::Y,
        2 => NodeAxis::Z,
        _ => return None,
    };
    {
        let p = &mut photons[begin..end];
        match axis {
            NodeAxis::X => {
                p.sort_by(|a, b| a.pos.x.partial_cmp(&b.pos.x).unwrap());
            }
            NodeAxis::Y => {
                p.sort_by(|a, b| a.pos.y.partial_cmp(&b.pos.y).unwrap());
            }
            NodeAxis::Z => {
                p.sort_by(|a, b| a.pos.z.partial_cmp(&b.pos.z).unwrap());
            }
        }
    }
    let median = (end - begin) / 2;

    let mut node = Box::new(TreeNode::new(photons[begin + median].clone(), axis));

    node.left_node = construct_kd_tree(photons, begin, begin + median, depth + 1);
    node.right_node = construct_kd_tree(photons, begin + median + 1, end, depth + 1);

    Some(node)
}
